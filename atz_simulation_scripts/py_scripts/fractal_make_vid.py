from transform import *

with open('three_to_one_300_1.02040816327.pdb', 'r') as myf:
    pdb_lines = myf.readlines()

#Need to block each model and keep appending them as new models
new_model = []
total = []
blank_addition = []
model_num = 0
for index, line in enumerate(pdb_lines):
    if "ENDMDL" in line:
        if len(blank_addition) < 55:
            blank_addition = [x[:21] + 'B' + x[22:] for x in blank_addition]
        else:
            blank_addition = [x[:21] + 'A' + x[22:] for x in blank_addition]
        total += blank_addition[:]
        blank_addition = []
        new_model.append(total[:])
    elif "MODEL" not in line:
        blank_addition.append(line)
    else:
        model_num += 1

#reconstruct the file in the right order
with open('fractal_growth.pdb', 'w') as newf:
    for index, model in enumerate(new_model, start=1):
        newf.write('MODEL  %s\n' % index)
        model = convert_to_pose_num(model)
        model = atom_renumber(model)
        newf.writelines(model)
        newf.write('ENDMDL\n')
        #convert_to_pose_num
        #atom_renum
