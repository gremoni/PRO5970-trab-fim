import pandas as pd
import io
import math
import gurobi as gp
from gurobipy import GRB 
from gurobipy import quicksum
import itertools
from itertools import permutations

url_hid = https://github.com/gremoni/PRO5970-trab-fim/blob/main/hidraulicas.csv
url_term = https://github.com/gremoni/PRO5970-trab-fim/blob/main/termicas.csv
url_cons = https://github.com/gremoni/PRO5970-trab-fim/blob/main/consumo.csv

csv_hid = pd.read_csv(url_hid, sep=';',  decimal=',', header=0, index_col=0)
csv_term = pd.read_csv(url_term, sep=';',header=0, index_col=0)
csv_cons = pd.read_csv(url_cons, sep=';', header=0, index_col=0)


# Transforma os CSVs em dicionários
dict_hid = csv_hid.to_dict()
dict_term = csv_term.to_dict()
dict_cons = csv_cons.to_dict()

# Cria listas auxiliares para as variáveis do modelo
list_hid = csv_hid.index.tolist()
list_term = csv_term.index.tolist()
list_cons = csv_cons.index.tolist()

# Separa as usinas por grupo proprietário 
gp_1 = ['h1', 'h2', 'h3', 'h4', 't1', 't2']
gp_2 = ['h5', 'h6', 'h7']
gp_3 = ['t3', 't4', 't5']

# Separa as usinas térmicas por grupo proprietário 
term_gp1 = [idx for idx in gp_1 if idx[0].lower() == 't'.lower()]
term_gp2 = [idx for idx in gp_2 if idx[0].lower() == 't'.lower()]
term_gp3 = [idx for idx in gp_3 if idx[0].lower() == 't'.lower()]

# Separa as usinas hidrelétricas por grupo proprietário 
hid_gp1 = [idx for idx in gp_1 if idx[0].lower() == 'h'.lower()]
hid_gp2 = [idx for idx in gp_2 if idx[0].lower() == 'h'.lower()]
hid_gp3 = [idx for idx in gp_3 if idx[0].lower() == 'h'.lower()]

# Separa as cascatas 
cascata_2 = ['h2', 'h3']
cascata_3 = ['h4', 'h7', 'h5']
cascata_4 = ['h6']


# Horizonte do problema (em meses)
T = 6

list_i = list()
for i in range(T):
    list_i.append(i+1)
    
beta = 50


################## DESPACHO POR CUSTO ##################

# Cria modelo no Gurobi
m = gp.Model()


# Cria as variáveis

GH = m.addVars(list_hid, list_i, vtype = GRB.CONTINUOUS, name='GH')
GT = m.addVars(list_term, list_i, vtype = GRB.CONTINUOUS, name='GT')
ear = m.addVars(list_hid, list_i, vtype = GRB.CONTINUOUS, name='ear')
vt = m.addVars(list_hid, list_i, vtype = GRB.CONTINUOUS, name='vt')
ve = m.addVars(list_hid, list_i, vtype = GRB.CONTINUOUS, name='ve')


# define a função objetivo
m.setObjective(
    quicksum(dict_term['ct'][t]*GT[t,i] for t in list_term for i in list_i)
    +quicksum(dict_hid['ch'][h]*GH[h,i] for h in list_hid for i in list_i),
    sense=gp.GRB.MINIMIZE)


    # define as restrições

#C1 = geração deve atender o consumo
c1 = m.addConstrs(gp.quicksum(GT[t,i] for t in list_term) 
                  + gp.quicksum(GH[h,i] for h in list_hid) 
                  == gp.quicksum(dict_cons['consumo'][c] for c in list_cons) for i in list_i)

#C2 = limites de geração
c2 = m.addConstrs(GH[h,i] <= dict_hid['ghmax'][h] for h in list_hid for i in list_i)
c3 = m.addConstrs(GT[t,i] <= dict_term['gtmax'][t] for t in list_term for i in list_i)

# C4
#c4 = m.addConstrs(ear[h,i] == ear[h,i-1]-((vt[h,i-1]+ve[h,i]-dict_hid['inc'][h]))/2.628e+6 for h in list_hid for i in list_i[1:])
c4 = m.addConstrs(ear[h,i] == ear[h,i-1]-((vt[h,i-1]-dict_hid['inc'][h]))/1e6 for h in list_hid for i in list_i[1:])
c5 = m.addConstrs(ear[h,1] == dict_hid['v1'][h] for h in list_hid)
c6 = m.addConstrs(ear[h,i] <= dict_hid['vmax'][h] for h in list_hid for i in list_i)
c7 = m.addConstrs(ear[h,i] >= dict_hid['vmin'][h] for h in list_hid for i in list_i)
c8 = m.addConstrs(ear[h,T] >=  dict_hid['vobj'][h] for h in list_hid)
c9 = m.addConstrs(GH[h,i] == dict_hid['rho'][h]*vt[h,i] for h in list_hid for i in list_i)


m.optimize()