# this script depends on sparta+
import os

batch_size = 128

protein_path = './pdbs'
output_path = './chemical_shifts'
current_path = os.getcwd()
sparta_path = 'ABSOLUTE/PATH/TO/SPARTA+'

# check if the output path exists
if not os.path.exists(output_path):
    os.makedirs(output_path)

protein_files = [os.path.join(protein_path, file)
                 for file in os.listdir(protein_path)
                 if file.endswith('.pdb')]

for iteration in range(0, len(protein_files), batch_size):
    # check current directory
    if not os.getcwd() == current_path:
        os.chdir(current_path)

    if iteration + batch_size < len(protein_files):
        batch_files = protein_files[iteration:iteration + batch_size]
    else:
        batch_files = protein_files[iteration:]

    # convert the file list into a string
    protein_files_str = ' '.join(batch_files)

    # generate command for SPARTA+
    with open('./chemical_shifts/cmd.sh', 'w') as f:
        f.write(f'{sparta_path} -in {protein_files_str}')

    # check if cmd.sh is execuble
    os.system('chmod +x ./chemical_shifts/cmd.sh')

    # get into the 'chemical_shifts' directory
    os.chdir('./chemical_shifts')

    # run the command
    os.system('./cmd.sh > log.txt')
    os.remove('log.txt')
