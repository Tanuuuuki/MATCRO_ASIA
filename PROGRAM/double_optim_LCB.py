import pandas as pd
import numpy as np
import GPyOpt
import GPy
import os
import sys
import time
import subprocess
import functools
print = functools.partial(print, flush=True)

random_seed = 100
np.random.seed(random_seed)

arg = sys.argv
country = arg[1]
region = arg[2]
year_st = int(arg[3])
year_ed = int(arg[4])
obj_path = arg[5].replace(' ', '')
iprm_path = arg[6].replace(' ', '')
yld_path = arg[7].replace(' ', '')
opt_path = arg[8].replace(' ', '')
step = arg[9].replace(' ', '')
season_flag = int(arg[10])
year_len = year_ed-year_st+1
season_list = ["1st", "2nd", "3rd"]

zrostep = "zero"
fststep = "first"
sndstep = "second"
trdstep = "third"

iniscript = f"./RUN/MATCRO_Asia ./SET/{country}/{region}/SET_{region}_INI.txt"
optscript = f"./RUN/MATCRO_Asia ./SET/{country}/{region}/SET_{region}_opt.txt"

# PRM ファイルの読み込み
if year_st < 1990:
    if step == fststep:
        iprm_file = f"{iprm_path}{year_st+10}-{year_ed+10}.txt"
    elif step == trdstep or step == sndstep:
        iprm_file = f"{iprm_path}{year_st}-{year_ed}.txt"
elif year_st >= 1990:
    if step == fststep:
        iprm_file = f"./SMPL/INI/{country}/SMPL_{region}.txt"
    elif step == trdstep or step == sndstep:
        iprm_file = f"{iprm_path}{year_st}-{year_ed}.txt"    
        
##### make phenology data #####
#igdh_file = f"./GDH/INI/GDH_{region}.txt"
#optgdh = pd.read_table(igdh_file, header=None)[0].values

for season_num in range(season_flag):
    season = season_list[season_num]
    
    if os.path.exists(iprm_file):
        optpar = pd.read_table(iprm_file, header=None)[0].values
    else:
        if year_ed >= 2015-20+1 and step == fststep:
            gdhm_file_ir = f"/home/tnakagawa/MATCRO_ASIA/INPUT/GDHm/{country}/{region}_ir_GDHm.csv"
            gdhm_ir = pd.read_csv(gdhm_file_ir)#, usecols=[1, 2, 3])
            optpar = [0.45, 1.0, 24.0, 37.0, 400]
        else :
            print("file does not exsit")
            print(year_st, step, iprm_file)
            os._exit(1)


with open(f"./SMPL/{country}/{region}/SMPL_{region}_INI.txt", 'w') as f:
    for value in optpar:
        f.write(str(value) + '\n')

##### PROGRAM #####
x = list(range(year_st, (year_ed + 1)))
nosui_file = f"./STATS/{country}_yield.csv"
nosui = pd.read_csv(nosui_file)
nosui_mod = nosui[nosui['region'] == region]
styear = min(nosui_mod['year'])
edyear = max(nosui_mod['year'])

obs_yld_set = [0] * (edyear-styear+1)
ir_fraction = [0] * (edyear-styear+1)
rf_fraction = [0] * (edyear-styear+1)

for i in range(len(nosui_mod['year'])):
    year = nosui_mod['year'][i]
    obs_yld_set[year-styear] = nosui_mod['yield_kg10a'][i] * 10

#print(obs_yld_set)

print(styear, year_st, year_ed)
obs_yld = obs_yld_set[(year_st-styear):(year_ed-styear+1)]
print(obs_yld)

###### parameter range #####
# 1 = SLN, 2 = HI, 3 = TCmin, 4 = THcrit, 5 = SLNY1

up1 = 0.8
up2 = 1.0
up3 = 30.0
up4 = 42.0
up5 = 500
up6 = 3000

lw1 = 0.0
lw2 = 0.5
lw3 = 15.0
lw4 = 30.0
lw5 = 250
lw6 = 800

#ir_file = f'./INPUT/IRR_RF/{country}/MOD/{region}_ir.csv'
#rf_file = f'./INPUT/IRR_RF/{country}/{region}_rf.csv'
#ir_frac = pd.read_csv(ir_file)[season]

##### Initial Yld Calculation #####
ini_yld = [0] * (year_ed-year_st+1)
for season_num in range(season_flag):
    season = season_list[season_num]
    iniscript_ir = f"./RUN/MATCRO_Asia ./SET/{country}/{region}/SET_{region}_ir{season}_INI.txt"
    iniscript_rf = f"./RUN/MATCRO_Asia ./SET/{country}/{region}/SET_{region}_rf{season}_INI.txt"
    
    subprocess.run(iniscript_ir, shell=True)
    subprocess.run(iniscript_rf, shell=True)

    print("end init")
    for irrf in ['ir', 'rf']:
        iyld_name = f'./YLD/{country}/{region}/YLD_{region}_{irrf}{season}_INI.txt'
        irrf_name = f'./INPUT/IRR_RF/{country}/MOD/{region}_{irrf}.csv'
        ini_yld_ori = pd.read_table(iyld_name, header=None)[0].values
        irrf_fraction = pd.read_csv(irrf_name)[season]
        ini_yld = ini_yld + ini_yld_ori * irrf_fraction[year_st-1896:year_ed-1896+1]
        #ini_prm = pd.read_table(iprm_name, header=None)[0].values

print(ini_yld)
print("wait for 5 seconds")
time.sleep(5)
#os._exit(1)

print('optimization start')

period = year_ed - year_st + 1

optscript_ir = f"./RUN/MATCRO_Asia ./SET/{country}/{region}/SET_{region}_ir{season}_opt.txt"
optscript_rf = f"./RUN/MATCRO_Asia ./SET/{country}/{region}/SET_{region}_rf{season}_opt.txt"

#oyld_name = f'./YLD/{country}/{region}/YLD_{region}_ir{season}_opt.txt'
#olai_name = f'./LAImx/{region}/LAImx_{region}_opt.txt'

# 最小化したい目的関数を定義
def objective_function(X):
    X = X[0]
    #SMPL = [round(X[0],1),round(X[1],2),round(X[2],1),round(X[3],3),round(X[4],1)]
    SMPL = [round(X[0],2),round(X[1],2),round(X[2],1),round(X[3],1),round(X[4],0)]#, round(X[5],3)]#, round(X[5],3)]
    print(SMPL)

    with open(f"./SMPL/{country}/{region}/SMPL_{region}_opt.txt", 'w') as f:
        for value in SMPL:
            f.write(str(value) + '\n')

    opt_yld = [0] * 20
    for season_num in range(season_flag):
        season = season_list[season_num]
        optscript_ir = f"./RUN/MATCRO_Asia ./SET/{country}/{region}/SET_{region}_ir{season}_opt.txt"
        optscript_rf = f"./RUN/MATCRO_Asia ./SET/{country}/{region}/SET_{region}_rf{season}_opt.txt"
    
        subprocess.run(optscript_ir, shell=True)
        subprocess.run(optscript_rf, shell=True)

        for irrf in ['ir', 'rf']:
            oyld_name = f'./YLD/{country}/{region}/YLD_{region}_{irrf}{season}_opt.txt'
            opt_yld_ori = pd.read_table(oyld_name, header=None)[0].values
            irrf_name = f'./INPUT/IRR_RF/{country}/MOD/{region}_{irrf}.csv'
            irrf_fraction = pd.read_csv(irrf_name)[season]
            opt_yld = opt_yld + opt_yld_ori*irrf_fraction[year_st-1896:year_ed-1896+1]

    yld_sum = np.sum(((obs_yld - opt_yld) ** 2)/len(obs_yld))
    yld_sq = np.sqrt(yld_sum)
    output = yld_sq

    print(output)
    return output
    

# 初期データ点を設定
if step == zrostep:
    space = [
        {'name': 'LEFY0',  'type': 'continuous', 'domain': (optpar[0], optpar[0])},
        {'name': 'PNCLY2', 'type': 'continuous', 'domain': (optpar[1], optpar[1])},
        {'name': 'TCmin',  'type': 'continuous', 'domain': (optpar[2], optpar[2])},
        {'name': 'THcrit', 'type': 'continuous', 'domain': (optpar[3], optpar[3])},
        {'name': 'SLWYA',  'type': 'continuous', 'domain': (optpar[4], optpar[4])},
        #{'name': 'GDHm',   'type': 'continuous', 'domain': (lw6, up6)}
    ]
elif step == fststep:
    space = [
        {'name': 'LEFY0',  'type': 'continuous', 'domain': (optpar[0], optpar[0])},
        {'name': 'PNCLY2', 'type': 'continuous', 'domain': (optpar[1], optpar[1])},
        {'name': 'TCmin',  'type': 'continuous', 'domain': (optpar[2], optpar[2])},
        {'name': 'THcrit', 'type': 'continuous', 'domain': (optpar[3], optpar[3])},
        {'name': 'SLWYA',  'type': 'continuous', 'domain': (lw5, up5)},
        #{'name': 'GDHm',   'type': 'continuous', 'domain': (optpar[5], optpar[5])}
    ]
elif step == sndstep:
    space = [
        {'name': 'LEFY0',  'type': 'continuous', 'domain': (lw1, up1)},
        {'name': 'PNCLY2', 'type': 'continuous', 'domain': (lw2, up2)},
        {'name': 'TCmin',  'type': 'continuous', 'domain': (optpar[2], optpar[2])},
        {'name': 'THcrit', 'type': 'continuous', 'domain': (optpar[3], optpar[3])},
        {'name': 'SLWYA',  'type': 'continuous', 'domain': (optpar[4], optpar[4])},
        #{'name': 'GDHm',   'type': 'continuous', 'domain': (optpar[5], optpar[5])}
    ]
elif step == trdstep:
    space = [
        {'name': 'LEFY0',  'type': 'continuous', 'domain': (optpar[0], optpar[0])},
        {'name': 'HI',     'type': 'continuous', 'domain': (optpar[1], optpar[1])},
        {'name': 'TCmin',  'type': 'continuous', 'domain': (lw3, up3)},
        {'name': 'THcrit', 'type': 'continuous', 'domain': (lw4, up4)},
        {'name': 'SLWYA',  'type': 'continuous', 'domain': (optpar[4], optpar[4])},
        #{'name': 'GDHm',   'type': 'continuous', 'domain': (lw6, up6)}
    ]


# 最適化を実行
initial_samples = np.array([[
    optpar[0],optpar[1],optpar[2],optpar[3],optpar[4]#,optpar[5]
]])

# Bayesian Optimization オブジェクトを作成
bo_optimizer = GPyOpt.methods.BayesianOptimization(
    f=objective_function, 
    domain=space, 
    X=initial_samples,
    acquisition_type='LCB'
    #constraints=constraints
    )
bo_optimizer.run_optimization(max_iter=10)

# 結果を表示
print('Optimization ended')
#print("random_seed:",random_seed)
print(bo_optimizer)
print("最適なパラメータ:", bo_optimizer.x_opt)
print("最小値:", bo_optimizer.fx_opt)

output = bo_optimizer.x_opt
output = [
    round(output[0],2),
    round(output[1],2),
    round(output[2],1),
    round(output[3],1),
    round(output[4],0),
    #round(output[5],3)
]

with open(f"./SMPL/{country}/{region}/SMPL_{region}_opt.txt", 'w') as f:
    for value in output :
        f.write(str(value) + '\n')

opt_yld = [0] * 20
opt_prm = [0] * 20
for season_num in range(season_flag):
    season = season_list[season_num]
    optscript_ir = f"./RUN/MATCRO_Asia ./SET/{country}/{region}/SET_{region}_ir{season}_opt.txt"
    optscript_rf = f"./RUN/MATCRO_Asia ./SET/{country}/{region}/SET_{region}_rf{season}_opt.txt"

    subprocess.run(optscript_ir, shell=True)
    subprocess.run(optscript_rf, shell=True)

    for irrf in ['ir', 'rf']:
        oyld_name = f'./YLD/{country}/{region}/YLD_{region}_{irrf}{season}_opt.txt'
        opt_yld_ori = pd.read_table(oyld_name, header=None)[0].values
        irrf_name = f'./INPUT/IRR_RF/{country}/MOD/{region}_{irrf}.csv'
        irrf_fraction = pd.read_csv(irrf_name)[season]
        opt_yld = opt_yld + opt_yld_ori*irrf_fraction[year_st-1896:year_ed-1896+1]        


ydat = pd.DataFrame({
    'year': x,
    'obs': obs_yld,
    'ini': ini_yld,
    'opt': opt_yld
})

optX = pd.DataFrame(
    bo_optimizer.X
)

optY = pd.DataFrame(
    bo_optimizer.Y
)

record_row = pd.DataFrame(range(1,(len(bo_optimizer.X)+1)))
aquisi_row = pd.DataFrame(range(1,(len(bo_optimizer.Y)+1)))

record = pd.concat([record_row, optX], axis=1)
aquisition = pd.concat([aquisi_row, optY], axis=1)

record.columns = 'year','LEFY0','PNCLY2','TCmin','THcrit','SLWYA'#,'GDHm'
aquisition.columns = 'year','aquisition'

ydat.to_csv(f"{yld_path}{year_st}-{year_ed}.csv", index=False)
record.to_csv(f"{obj_path}/OPT_TRACE_{year_st}-{year_ed}.csv", index=False)
aquisition.to_csv(f"{obj_path}/OPT_ACQ_{year_st}-{year_ed}.csv", index=False)

with open(f"{opt_path}{year_st}-{year_ed}.txt", 'w') as f:
    for value in output :
        f.write(str(value) + '\n')

sys.exit()
