import os

import numpy as np

target_shift = 'N'
target_method = 'chemical'

protein_path = f'./{target_method}_shifts'
exp_path = './cs_ref'
output_path = './idpfold_cs_output'
protein_name_list = [i.split('.')[0] for i in os.listdir(f'{exp_path}/{target_shift}')]

if not os.path.exists(output_path):
    os.makedirs(output_path)

average_cs = []
exp_cs = []
for protein_name in protein_name_list:

    # calculate experimental chemical shifts
    exp_cs_index = []
    exp_cs_value = []
    with open(f'{exp_path}/{target_shift}/{protein_name}.dat', 'r') as f:

        for lines in f.readlines():
            if not lines.startswith('#'):
                exp_cs_index.append(int(lines.split()[0][1:]))
                exp_cs_value.append(float(lines.split()[1]))

    current_cs = []

    protein_files = [os.path.join(protein_path, file)
                     for file in os.listdir(protein_path)
                     if file.endswith('_pred.tab') and file.startswith(protein_name)]

    for files in protein_files:

        file_index = 0
        current_conf = []
        current_exp_cs = []
        detected = False

        if target_shift == 'N':
            current_conf.append(0.)
            current_exp_cs.append(0.)

        with open(files, 'r') as f:
            file_content = f.readlines()
        for lines in file_content:
            if detected is False:
                if lines.startswith('FORMAT'):
                    detected = True
                continue

            if len(lines.split()) >= 6:
                if target_shift == 'CB':
                    if lines.split()[2] == 'CA':
                        current_conf.append(0.)

                        if file_index == 0:
                            current_exp_cs.append(0.)
                    if lines.split()[2] == 'CB':
                        current_conf[-1] += float(lines.split()[3])

                        if file_index == 0:
                            current_exp_cs[-1] += float(lines.split()[5])
                elif target_shift == 'HN':
                    if lines.split()[2] == 'CA':
                        current_conf.append(0.)

                        if file_index == 0:
                            current_exp_cs.append(0.)
                    if lines.split()[2] == 'HN':
                        current_conf[-1] += float(lines.split()[3])

                        if file_index == 0:
                            current_exp_cs[-1] += float(lines.split()[5])
                elif target_shift == 'HA':
                    if lines.split()[2] == 'HA' or lines.split()[2] == 'HA2':
                        current_conf.append(float(lines.split()[3]))

                        if file_index == 0:
                            current_exp_cs.append(float(lines.split()[5]))
                elif target_shift == 'N':
                    if lines.split()[2] == 'N':
                        current_conf.append(float(lines.split()[3]))

                        if file_index == 0:
                            current_exp_cs.append(float(lines.split()[5]))

                    if lines.split()[2] == 'CA' and lines.split()[1] == 'P':
                        current_conf.append(0.)

                        if file_index == 0:
                            current_exp_cs.append(0.)
                elif target_shift == 'CA':
                    if lines.split()[2] == 'CA':
                        current_conf.append(float(lines.split()[3]))

                        if file_index == 0:
                            current_exp_cs.append(float(lines.split()[5]))

        file_index += 1

        helix_prop = np.zeros(len(current_conf))
        for i in range(len(current_conf)):
            if current_conf[i] > 0.5:
                helix_prop[i] = 1

        current_cs.append(current_conf)

    current_cs = np.array(current_cs)
    current_avg_cs = np.nanmean(current_cs, axis=0)

    for index, value in zip(exp_cs_index, exp_cs_value):
        if protein_name == 'Sic1':
            current_exp_cs[index - 3] = value - current_exp_cs[index - 3]
        else:
            current_exp_cs[index - 1] = value - current_exp_cs[index - 1]

    # change all index not in exp_cs_index to np.nan
    for i in range(len(current_avg_cs)):
        if i not in exp_cs_index:
            current_exp_cs[i - 1] = np.nan

    average_cs.append(current_avg_cs)
    exp_cs.append(current_exp_cs)

# make a dictionary
cs_dict = {protein_name: cs for protein_name, cs in zip(protein_name_list, average_cs)}
np.save(f'{output_path}/{target_method}_dict_{target_shift}.npy', cs_dict)

exp_dict = {protein_name: cs for protein_name, cs in zip(protein_name_list, exp_cs)}
np.save(f'{output_path}/exp_dict_{target_shift}.npy', exp_dict)

# calculate rmsd
rmsd = []
for protein_name in protein_name_list:
    rmsd.append(np.sqrt(np.nanmean((cs_dict[protein_name] - exp_dict[protein_name]) ** 2)))

print(f'Mean RMSD: {np.mean(rmsd):.3f}')
