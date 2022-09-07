def oti_gp1():
    global COFT
    global COFH
    global QOH
    global QOT
    global GH 
    global GT 
    global ear
    global vt 
    global ve 
    global spot

    #OTIMIZANDO PARA GRUPO 1
    
    ################## DESPACHO POR OFERTA ##################
    # Cria modelo no Gurobi
    m = gp.Model()
    
    # Cria as variáveis
    GH = m.addVars(list_hid, list_i, vtype = GRB.CONTINUOUS, name='GH')
    GT = m.addVars(list_term, list_i, vtype = GRB.CONTINUOUS, name='GT')
    ear = m.addVars(list_hid, list_i, vtype = GRB.CONTINUOUS, name='ear')
    vt = m.addVars(list_hid, list_i, vtype = GRB.CONTINUOUS, name='vt')
    ve = m.addVars(list_hid, list_i, vtype = GRB.CONTINUOUS, name='ve')
    
    spot = m.addVars(list_i, vtype = GRB.CONTINUOUS, name='spot')
    DFGH = m.addVars(list_hid, list_i, vtype = GRB.CONTINUOUS, name='DFGH')
    QOH = m.addVars(list_hid, list_i, vtype = GRB.CONTINUOUS, name='QOH')
    QOT = m.addVars(list_term, list_i, vtype = GRB.CONTINUOUS, name='QOT')
    COFT = m.addVars(list_term, list_i, vtype = GRB.CONTINUOUS, name='COFT')
    COFH = m.addVars(list_hid, list_i, vtype = GRB.CONTINUOUS, name='COFH')
    delta_t_min = m.addVars(list_term, list_i, vtype = GRB.CONTINUOUS, name='delta_t_min')
    delta_t_max = m.addVars(list_term, list_i, vtype = GRB.CONTINUOUS, name='delta_t_max')
    delta_h_min = m.addVars(list_hid, list_i, vtype = GRB.CONTINUOUS, name='delta_h_min')
    delta_h_max = m.addVars(list_hid, list_i, vtype = GRB.CONTINUOUS, name='delta_h_max')
    
    
    # define a função objetivo
    m.setObjective(
        quicksum((spot[i]-dict_term['ct'][t])*GT[t,i] for t in term_gp1 for i in list_i)
        +quicksum((spot[i]-dict_hid['ch'][h])*GH[h,i] for h in hid_gp1 for i in list_i)
        -beta*quicksum(DFGH[h,i] for h in hid_gp1 for i in list_i),
        sense=gp.GRB.MAXIMIZE)

    
    

    # define as restrições
    #primeiro nível
    c12_1 = m.addConstrs(QOH[h,i] <= dict_hid['ghmax'][h] for h in hid_gp1 for i in list_i)
    c12_2 = m.addConstrs(QOH[h,i] >= 0 for h in hid_gp1 for i in list_i)
    c13_1 = m.addConstrs(QOT[t,i] <= dict_term['gtmax'][t] for t in term_gp1 for i in list_i)
    c13_2 = m.addConstrs(QOT[t,i] >= 0 for t in term_gp1 for i in list_i)
    c14 = m.addConstrs(ear[h,i] == 
                       ear[h,i-1]-((vt[h,i-1]+ve[h,i]-dict_hid['inc'][h]))/1e6 
                       for h in hid_gp1 for i in list_i[1:])
    c15 = m.addConstrs(GH[h,i] - dict_hid['rho'][h]*vt[h,i] <= DFGH[h,i]  for h in hid_gp1 for i in list_i)
    c16 = m.addConstrs(-(GH[h,i] - dict_hid['rho'][h]*vt[h,i]) <= DFGH[h,i]  for h in hid_gp1 for i in list_i)
    c17 = m.addConstrs(DFGH[h,i] >= 0 for h in hid_gp1 for i in list_i)
    c18_1 = m.addConstrs(ear[h,i] <= dict_hid['vmax'][h] for h in hid_gp1 for i in list_i)
    c18_2 = m.addConstrs(ear[h,i] >= dict_hid['vmin'][h] for h in hid_gp1 for i in list_i)
    c18_3 = m.addConstrs(ear[h,T] >=  dict_hid['vobj'][h] for h in hid_gp1)

    #segundo nível
    c23 = m.addConstrs(gp.quicksum(GT[t,i] for t in list_term) 
                       + gp.quicksum(GH[h,i] for h in list_hid) 
                       == gp.quicksum(dict_cons['consumo'][c] for c in list_cons) for i in list_i)
    c24_1 = m.addConstrs(GT[t,i] <= QOT[t,i] for t in list_term for i in list_i)
    c24_2 = m.addConstrs(GT[t,i] >= 0 for t in list_term for i in list_i)
    c25_1 = m.addConstrs(GH[h,i] <= QOH[h,i] for h in list_hid for i in list_i)
    c25_2 = m.addConstrs(GH[h,i] >= 0 for h in list_hid for i in list_i)
    c26 = m.addConstrs(spot[i] - delta_t_min[t,i] + delta_t_max[t,i] <= COFT[t,i] for t in list_term for i in list_i)
    c27 = m.addConstrs(spot[i] - delta_h_min[h,i] + delta_h_max[h,i] <= COFH[h,i] for h in list_hid for i in list_i)
    c28_1 = m.addConstrs(delta_t_min[t,i] <= 0 for t in list_term for i in list_i)
    c28_2 = m.addConstrs(delta_t_max[t,i] <= 0 for t in list_term for i in list_i)
    c28_3 = m.addConstrs(delta_h_min[h,i] <= 0 for h in list_hid for i in list_i)
    c28_4 = m.addConstrs(delta_h_max[h,i] <= 0 for h in list_hid for i in list_i)
    c5 = m.addConstrs(ear[h,1] == dict_hid['v1'][h] for h in list_hid)

    
    
    #c29 = m.addConstrs(GT[t,i]*(spot[i]-delta_t_min[t,i]+delta_t_max[t,i]-COFT[t,i])==0 for t in list_term for i in list_i)
    #c30 = m.addConstrs(GH[h,i]*(spot[i]-delta_h_min[h,i]+delta_h_max[h,i]-COFH[h,i])==0 for h in list_hid for i in list_i)
    #c31 = m.addConstrs(spot[i]*(gp.quicksum(GH[h,i] for h in list_hid) + gp.quicksum(GT[t,i] for t in list_term) - gp.quicksum(dict_cons['consumo'][c] for c in list_cons))==0 for i in list_i)
    #c32 = m.addConstrs(delta_t_min[t,i]*GT[t,i] == 0 for t in list_term for i in list_i)
    #c33 = m.addConstrs(delta_h_min[h,i]*GH[h,i] == 0 for h in list_hid for i in list_i)
    #c34 = m.addConstrs(delta_t_max[t,i]*(QOT[t,i]-GT[t,i]) == 0 for t in list_term for i in list_i)
    #c35 = m.addConstrs(delta_h_max[h,i]*(QOH[h,i]-GH[h,i]) == 0 for h in list_hid for i in list_i)


    c35 = m.addConstr(gp.quicksum(COFT[t,i]*GT[t,i] for t in list_term for i in list_i)
                       + gp.quicksum(COFH[h,i]*GH[h,i] for h in list_hid for i in list_i)
                      -gp.quicksum(dict_cons['consumo'][c]*spot[i] for c in list_cons for i in list_i)
                      -gp.quicksum(QOT[t,i]*delta_t_max[t,i] for t in list_term for i in list_i)
                      -gp.quicksum(QOH[h,i]*delta_h_max[h,i] for h in list_hid for i in list_i)
                      == 0)

    
    #conhecimento das ofertas dos demais players
    c_teste1 = m.addConstrs(QOT[t,i] == QOT_gp2_ant[t,str(i)] for t in term_gp2 for i in list_i)
    c_teste12 = m.addConstrs(QOT[t,i] == QOT_gp3_ant[t,str(i)] for t in term_gp3 for i in list_i)
    c_teste2 = m.addConstrs(QOH[h,i] == QOH_gp3_ant[h,str(i)] for h in hid_gp3 for i in list_i)
    c_teste21 = m.addConstrs(QOH[h,i] == QOH_gp2_ant[h,str(i)] for h in hid_gp2 for i in list_i)
    c_teste3 = m.addConstrs(COFT[t,i] == COFT_gp3_ant[t,str(i)] for t in term_gp3 for i in list_i)
    c_teste31 = m.addConstrs(COFT[t,i] == COFT_gp2_ant[t,str(i)] for t in term_gp2 for i in list_i)
    c_teste4 = m.addConstrs(COFH[h,i] == COFH_gp3_ant[h,str(i)] for h in hid_gp3 for i in list_i)
    c_teste41 = m.addConstrs(COFH[h,i] == COFH_gp2_ant[h,str(i)] for h in hid_gp2 for i in list_i)
    
    #GurobiError: Q matrix is not positive semi-definite (PSD). Set NonConvex parameter to 2 to solve model.
    m.params.NonConvex = 2
    m.params.TimeLimit = 180
    
    m.optimize()
    
    #CUSTOS DAS OFERTAS TÉRMICAS ANTERIORES
    global COFT_gp1_ant
    global  QOT_gp1_ant
    global COFH_gp1_ant
    global  QOH_gp1_ant
    
    COFT_gp1_ant = {}
    
    for t in term_gp1:
        for i in (range(1,T+1)):
            COFT_gp1_ant["{}".format(t),  '{}'.format(i)] = COFT.select(t, i)[0].X
    
            
    # QUANTIDADES DAS OFERTAS TÉRMICAS ANTERIORES
    QOT_gp1_ant = {}
    
    for t in term_gp1:
        for i in range(1,T+1):
            QOT_gp1_ant["{}".format(t), "{}".format(i)] = QOT.select(t, i)[0].X
            
    #CUSTOS DAS OFERTAS HIDRELÉTRICAS ANTERIORES
    COFH_gp1_ant = {}
    
    for h in hid_gp1:
        for i in (range(1,T+1)):
            COFH_gp1_ant["{}".format(h),  '{}'.format(i)] = COFH.select(h, i)[0].X
    
            
    # QUANTIDADES DAS OFERTAS HIDRELÉTRICAS ANTERIORES
    QOH_gp1_ant = {}
    
    for h in hid_gp1:
        for i in range(1,T+1):
            QOH_gp1_ant["{}".format(h), "{}".format(i)] = QOH.select(h, i)[0].X
              
                        
    global spot_var
    spot_var = {}
    for i in range(1,T+1):
        spot_var["{}".format(i)] = spot.select(i)[0].X 
        
    global L_g1
    L_g1 = m.objVal
       
    
    global GT_total
    global GH_total
    global GH_total_gp1
    global GH_total_gp2
    global GT_total_gp1
    global GT_total_gp3

    GT_total = {}
    GH_total = {}
    
    GH_total_gp1=0
    GH_total_gp2=0
    GT_total_gp1=0
    GT_total_gp3=0

    
    for t in list_term:
        for i in range(1,T+1):
            GT_total["{}".format(t), "{}".format(i)] = GT.select(t, i)[0].X 
            
    for h in list_hid:
        for i in range(1,T+1):
            GH_total["{}".format(h), "{}".format(i)] = GH.select(h, i)[0].X 
            
    global ear_oferta 
    ear_oferta = {}
    
    for h in list_hid:
        for i in range(1,T+1):
            ear_oferta["{}".format(h), "{}".format(i)] = ear.select(h, i)[0].X