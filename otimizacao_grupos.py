i = 0

while (
    (COFH_gp1_ant != COFH_gp1_atu) or
    (COFH_gp2_ant != COFH_gp2_atu) or
    (COFT_gp1_ant != COFT_gp1_atu) or
    (COFT_gp2_ant != COFT_gp1_atu)):
       
       ofertas_inic()
       oti_gp1()
       COFH_gp1_atu = COFH_gp1_ant
       COFT_gp1_atu = COFT_gp1_ant

       oti_gp2()
       COFH_gp2_atu = COFH_gp2_ant

       oti_gp3()
        COFT_gp3_atu = COFT_gp3_ant
    
         i = i+1

    