#!/usr/bin/env python
import math, numpy, csv, sys, os, time, ast

'''
##### CORE ###############################################
## 
##  This file contains a suite of core functions that
##  are used to manipulate pdb files.
##  Author: William Hansen Email: wah49@scarletmail.rutgers.edu
##
##  
##  INDEX
##      1.Read_Write_Copy
##          copy_transform_object()
##          read_file()
##          read_write_file()
##          write_file()
##          read_dictionary_file()
##          write_dictionary_file()
##          pdb_to_xyz()
##          pose_to_xyz()
##          create_arbPDB_from_xyz()
##      2.PDB_slicing
##          clean_pdb()
##          atom_coordinate()
##          list_of_resAtoms()
##          list_of_residue_names()
##          list_of_chains()
##          convert_resname()
##          replace_coordLine()
##          grab_coords()
##          symmetric_chain_separate()
##          gen_residue_list()
##          block_pdb_by_res()
##          unblock_pdb_by_res()
##          convert_to_pose_num()
##          identify_pose_number()
##          convert_grid_string_to_list()
##          find_longest_grid_diameter()
##          atom_renumber()
##          combine_all_chains()
##          convert_to_coarse_model()
##      3.PDB_comparisons   
##          geometric_peptide_angle_check()
##          clash_check()  **LEGACY**
##      4.Vector_math
##          get_angle()
##          vector_angle()
##          get_mag()
##          get_mag_difference()
##          get_vector_difference()
##          get_dihedral()
##          normalize()
##          cross()
##          dot()
##          q_mult()
##          q_conjugate()
##          q_xyz_mult()
##          axisangle_to_q()
##      5.Transform_classes 
##          Transform()
##              get_xyz_info()
##              get_xyz_from_name_list()
##              update_coordinates_from_truncated_obj()
##              merge_secondary_object()
##              get_geo_center()
##              get_symmetric_axes()
##              get_xyz_by_chain()
##              get_xyz_by_atomtype()
##              get_xyz_by_restype()
##              get_xyz_by_resnum()
##              get_xyz_by_grid()
##              rename_pdblines_from_xyz_names()
##              generate_sub_symm_objs()
##              calculate_symmetric_transform()
##              rotate_object()
##              rotate_about_axis()
##              translate()
##              translate_along_axis()
##              generate_axes()
##              major_symmetric_axis_align()
##              align_object_to_vector()
##              get_transform_operations()
##              apply_transform_operations()
##              apply_symmetric_transform()
##              get_atom_grid()
##              clash_check()
##              return_to_pdb()
##              write_to_file()
##              append_to_file()
##          TransformCopy()
##              assign()
##          TransformVector()
##              rotate_vector()
'''

#### 1.Read_Write_Copy ###################################
##
##  Author: William Hansen

def copy_transform_object(obj):
    new_grid_name_dict = {}
    for k,v in obj.grid_name_dict.iteritems():
        new_grid_name_dict[k] = v[:]
    
    matrix_transform = {K:V[:] for K,V in obj.matrix_transform.iteritems()}
    sub_symm_objs = {k:copy_transform_object(v) for k,v in obj.sub_symm_objs.iteritems()}
    sub_symm_centers = {k:v[:] for k,v in obj.sub_symm_centers.iteritems()}

    info = {'pdb' : obj.pdb[:],\
            'xyz_total' : obj.xyz_total[:],\
            'xyz_names' : obj.xyz_names[:],\
            'geo_center_mags' : obj.geo_center_mags[:],\
            'grid_name_dict' : new_grid_name_dict,\
            'symm_type' : obj.symm_type,\
            'major_symm_axis' : obj.major_symm_axis,\
            'subunits' : obj.subunits,\
            'matrix_transform' : matrix_transform,\
            'chain_names' : obj.chain_names[:],\
            'sub_symm_objs' : sub_symm_objs,\
            'sub_symm_centers' : sub_symm_centers}
    new_obj = TransformCopy(info)
    return new_obj

def read_file(file_path):
    with open(file_path, 'r') as f:
        my_file = f.readlines()
    return my_file

def read_write_file(file_path):
    with open(file_path, 'r+') as f:
        my_file = f.readlines()
    return my_file

def write_file(input_pdb, path='unamed_pdb.pdb', chain=''):
        
    if not chain:
        with open(path, 'w') as myfile:
            myfile.writelines(input_pdb)
    else:
        pdb = clean_pdb(input_pdb)
        if type(chain) == list or type(chain) == tuple:
            pdb = [x for x in pdb if str(x[21]) in chain]
        elif type(chain) == str:
            pdb = [x for x in pdb if str(x[21]) == chain]
        with open(path, 'w') as myfile:
            myfile.writelines(pdb)

def read_dictionary_file(my_file):
    with open(my_file, 'r') as f:
        s = f.read()
        my_dict = ast.literal_eval(s)
    return my_dict

def write_dictionary_file(my_dict, my_file_name):
    with open(my_file_name, 'a') as mf:
        mf.write(str(my_dict))

def pdb_to_xyz(file_input): #was xyz_matrix_gen
    xyz = []
    names = []
    ##print("")
    for line in file_input:
        if line.startswith(("ATOM","HETATM")):
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            try:
                names.append(str(int(line[22:30])).upper()+ '_' +\
                                line[21].upper() + '_' +\
                                line[17:20].strip().upper() + '_' +\
                                line[11:17].strip().upper())
                #if line[12:16].strip() == "CA" or line[17:21].strip() == "****":
                xyz.append([x, y, z])
            except ValueError:
                continue
    
    xyz_coord = numpy.array(xyz[:], float)
    xyz_center = numpy.mean(xyz_coord, 0)
    #center each file with geometric center
    delta = list(xyz_coord - xyz_center)
    mags = [ sum(x**2)**0.5 for x in delta ]
    
    return xyz, names, mags

'''

def pose_to_xyz(pose):
    #rosetta.init()
    xyz = []
    total_res = pose.n_residue()
    for res in xrange(total_res):
        for atom in xrange(len(pose.residue(res+1))):
            xyz.append( [pose.residue(res+1)[atom][0], pose.residue(res+1)[atom][1], pose.residue(res+1)[atom][2]] )
    
    return xyz

def xyz_to_pose
    unfinished part of code, will get around to it

'''


def create_arbPDB_from_xyz(coords, file_name='', chain_id='A'): #was rep_coords
    
    preffix = "HETATM    1  X   LIG %s   1    " % chain_id
    suffix  = "  1.00  0.00           Z\n"
    out_pdb = []
    for i, coord in enumerate(coords):
        if abs(coord[0]) >= 100. and abs(coord[0]) < 1000. and coord[0] < 0.:
            x_comp = "%.2f" % coord[0]
        elif abs(coord[0]) >= 1000. and abs(coord[0]) < 10000.:
            x_comp = "%.2f" % coord[0]
            if coord[0] < 0.:
                x_comp = "%.1f" % coord[0]
        elif abs(coord[0]) >= 10000. and abs(coord[0]) < 100000.:
            x_comp = "%.1f" % coord[0]
            if coord[0] < 0.:
                x_comp = "%.0f" % coord[0]
        elif abs(coord[0]) >= 100000. and abs(coord[0]) < 1000000.:
            x_comp = int(coord[0])
        else:
            x_comp = "%.3f" % coord[0]

        if abs(coord[1]) >= 100. and abs(coord[1]) < 1000. and coord[1] < 0.:
            y_comp = "%.2f" % coord[1]
        elif abs(coord[1]) >= 1000. and abs(coord[1]) < 10000.:
            y_comp = "%.2f" % coord[1]
            if coord[1] < 0.:
                y_comp = "%.1f" % coord[1]
        elif abs(coord[1]) >= 10000. and abs(coord[1]) < 100000.:
            y_comp = "%.1f" % coord[1]
            if coord[1] < 0.:
                y_comp = "%.0f" % coord[1]
        elif abs(coord[1]) >= 100000. and abs(coord[1]) < 1000000.:
            y_comp = int(coord[1])
        else:
            y_comp = "%.3f" % coord[1]

        if abs(coord[2]) >= 100. and abs(coord[2]) < 1000. and coord[2] < 0.:
            z_comp = "%.2f" % coord[2]
        elif abs(coord[2]) >= 1000. and abs(coord[2]) < 10000.:
            z_comp = "%.2f" % coord[2]
            if coord[2] < 0.:
                z_comp = "%.1f" % coord[2]
        elif abs(coord[2]) >= 10000. and abs(coord[2]) < 100000.:
            z_comp = "%.1f" % coord[2]
            if coord[2] < 0.:
                z_comp = "%.0f" % coord[2]
        elif abs(coord[2]) >= 100000. and abs(coord[2]) < 1000000.:
            z_comp = int(coord[2])
        else:
            z_comp = "%.3f" % coord[2]


        line_out = ''.join([preffix, ('{0:>8}').format(x_comp), \
                                     ('{0:>8}').format(y_comp), \
                                     ('{0:>8}').format(z_comp), \
                            suffix])
        out_pdb.append(line_out)
    
    if file_name:
        write_file(out_pdb, file_name)
    
    return out_pdb


#### 2.PDB_slicing ###################################
##
##  Author: William Hansen and Elliott Dolan

def clean_pdb(pdb, add_headers=[]):
    keep_headers = ['ATOM  ','HETATM']
    if add_headers:
        if type(add_headers) == str:
            keep_headers.append(add_headers)
        else:
            for header in add_headers:
                keep_headers.append(header)
    return [x for x in pdb if x[0:6] in keep_headers]

#if no chain is given it assumes you mean the first name/residue number found
def atom_coordinate(pdb, name, residue=1, chain=[]):
    pdb = clean_pdb(pdb)
    for line in pdb:
        if int(line[22:26]) == int(residue):
            ##print 'hi'
            if type(name) == str:
                name = [name]
            if line[12:16].strip(' ') in name:
                ##print 'low'
                if not chain:
                    return [ float(line[30:38]),
                             float(line[38:46]),
                             float(line[46:54]), 
                             str(line[21]) ] 
                else:
                    if str(line[21]) == chain:
                        return [ float(line[30:38]),
                                 float(line[38:46]),
                                 float(line[46:54]),
                                 str(line[21]) ]

#if chain not specified then it will include all the atoms of all chains
def list_of_resAtoms(pdb, residue, chain=[]):
    pdb = clean_pdb(pdb)
    atom_names = []
    for line in pdb:
        if int(line[22:26]) == int(residue):
            ##print 'hi'
            if not chain:
                if line[11:17].strip(' ') not in atom_names:
                    atom_names.append(line[11:17].strip(' '))
            else:
                if str(line[21]) == chain:
                    atom_names.append(line[11:17].strip(' '))
    return atom_names

def list_of_residue_names(pdb):
    residue_names = []
    pdb = clean_pdb(pdb)
    pdb = block_pdb_by_res(pdb)
    for block in pdb:
        residue_names.append(block[0][17:20])
    return residue_names

def list_of_chains(pdb):
    chains = set()
    pdb = clean_pdb(pdb)
    for line in pdb:
        chains.add(line[21])
    chains = list(chains)
    chains.sort()
    return chains[:]

def convert_resname(res):
    
    res1_dictionary = {"A" : "ALA",\
                     "R" : "ARG",\
                     "N" : "ASN",\
                     "D" : "ASP",\
                     "C" : "CYS",\
                     "Q" : "GLN",\
                     "E" : "GLU",\
                     "G" : "GLY",\
                     "H" : "HIS",\
                     "I" : "ILE",\
                     "L" : "LEU",\
                     "K" : "LYS",\
                     "M" : "MET",\
                     "F" : "PHE",\
                     "P" : "PRO",\
                     "S" : "SER",\
                     "T" : "THR",\
                     "W" : "TRP",\
                     "Y" : "TYR",\
                     "V" : "VAL"}
    res3_dictionary = {v: k for k, v in res1_dictionary.iteritems()}
    if res in res1_dictionary.keys():
        return res1_dictionary[res]
    if res in res3_dictionary.keys():
        return res3_dictionary[res]
 
def replace_coordLine(line, coord):
    newline = ''.join([line[0:30], ('{0:>8}').format("%.3f" % coord[0]), \
                                   ('{0:>8}').format("%.3f" % coord[1]), \
                                   ('{0:>8}').format("%.3f" % coord[2]), line[54:]])
    return newline

def grab_coords(line):
    return [ float(line[30:38]),
             float(line[38:46]),
             float(line[46:54]) ] 

def symmetric_chain_separate(pdb, HETATM=False):
    chain_types = set()
    if HETATM == False:
        master = []
        for line in pdb:
            if line[0:6] in ["ATOM  "]:
                master.append(line)
                chain_types.add(line[21])
    else:
        master = []
        for line in pdb:
            if line[0:6] in ["ATOM  ", "HETATM"]:
                master.append(line)
                chain_types.add(line[21])
    subunits = len(chain_types)
    chains = zip(*(iter(master),)*(len(master)/subunits))
    return chains

def gen_residue_list(pdb):
    pdb = clean_pdb(pdb)
    res_nums = []
    for line in pdb:
        if int(line[22:26]) not in res_nums:
            res_nums.append(int(line[22:26]))
    return res_nums

def reletter_chain(pdb):
    pdb_block = []
    indiv_block = []
    for line in pdb:
        if line == pdb[0]:
            last_chain = line[21]
        current_chain = line[21]
        if current_chain == last_chain:
            indiv_block.append(line)
        else:
            pdb_block.append(indiv_block)
            indiv_block = []
            indiv_block.append(line)
            last_chain = line[21]
    pdb_block.append(indiv_block)
    letters = ['A','B','C','D','E','F','G','H',\
                'I','J','K','L','M','N','O','P',\
                'Q','R','S','T','U','V','W','X','Y','Z']
    relettered_pdb = []
    for index, chain in enumerate(pdb_block):
        for line in chain:
            new_line = line[0:21] + letters[index] + line[22:]
            relettered_pdb.append(new_line)
    return relettered_pdb

def block_pdb_by_res(pdb):
    pdb_block = []
    indiv_block = []
    for line in pdb:
        if line == pdb[0]:
            last_res = int(line[22:26])
        current_res = int(line[22:26])
        if current_res == last_res:
            indiv_block.append(line)
        else:
            pdb_block.append(indiv_block)
            indiv_block = []
            indiv_block.append(line)
            last_res = int(line[22:26])
    #append last block
    pdb_block.append(indiv_block)
    return pdb_block

def unblock_pdb_by_res(pdb_block):
    pdb = []
    for block in pdb_block:
        for line in block:
            pdb.append(line)
    return pdb

#must be a monomer to do this
def convert_to_pose_num(pdb):
    #current_res_nums = gen_residue_list(pdb)
    header = [x for x in pdb if x.startswith("REMARK")]
    pdb = clean_pdb(pdb)
    res_separated_blocks = block_pdb_by_res(pdb)
    fixed_blocks = []
    for i_block, block in enumerate(res_separated_blocks):
        new_block = []
        for i_line, line in enumerate(block):
            #force the line at residue position to change to the block index +1
            newline = ''.join([res_separated_blocks[i_block][i_line][:22], \
                      ('{0:>4}').format( "%i" % (i_block + 1) ), \
                      res_separated_blocks[i_block][i_line][26:]])
            new_block.append(newline)
        fixed_blocks.append(new_block)

    pdb = unblock_pdb_by_res(fixed_blocks)
    if not header:
        return pdb
    else:
        return header + pdb

def identify_pose_number(pdb, residue_num, chain='A'):
    pdb_obj = Transform(pdb)
    pdb_pose = convert_to_pose_num(pdb)
    pdb_pose_obj = Transform(pdb_pose)
    for i, name in enumerate(pdb_obj.xyz_names):
        if name.split('_')[1] == chain:
            if name.split('_')[0] == residue_num:
                pose_num = pdb_pose_obj.get_xyz_info(i)[1].split('_')[0]
                return pose_num

def convert_grid_string_to_list(grid):
    grid_list = [ast.literal_eval(grid_string) for grid_string in grid]
    return grid_list 

def find_longest_grid_diameter(grid):
    #grid_list = convert_grid_string_to_list(grid)
    mags = []
    for xyz1 in grid:
        for xyz2 in grid:
            mag = get_mag_difference(ast.literal_eval(xyz1), ast.literal_eval(xyz2))
            mags.append(mag)
    return max(mags)

def atom_renumber(pdb):
    pdb = clean_pdb(pdb)
    new_pdb = []
    atom_count = 0
    for lines in pdb:
        atom_count += 1
        new_pdb.append(lines[:6]+str(atom_count).rjust(5,' ') + lines[11:])
    return new_pdb

def combine_all_chains(pdb, chain_name='A'):
    header = [x for x in pdb if x.startswith(("REMARK"))]
    pdb = clean_pdb(pdb)
    out_pdb = []
    for line in pdb:
        newline = ''.join([line[:21], chain_name, line[22:]])
        out_pdb.append(newline)

    if not header:
        return out_pdb
    else:
        return header + out_pdb

#Convert to coarse takes a input number of points to bin the protein CA's and create a coarse model
## This currently does not work with asymmetric pdbs

def convert_to_coarse_model(pdb, cluster_num=1, by_chain=True):
    if type(cluster_num) == float:
        cluster_num = int(cluster_num)
    if cluster_num < 1:
        cluster_num = 1
    
    pdb_obj = Transform(pdb)
    orig_center = pdb_obj.get_geo_center()[:]
    pdb_obj.translate([0.,0.,0.], update=True)
    pdb_obj.get_xyz_by_atomtype('CA', update=True)
    pdb_obj.calculate_symmetric_transform()
    if by_chain == False:
        xyz_total_copy = pdb_obj.xyz_total[:]
    else:
        chain_list = list_of_chains(pdb)
        pdb_obj_main_chain = copy_transform_object(pdb_obj)
        xyz_total_copy = pdb_obj_main_chain.get_xyz_by_chain(chain_list[0])[:]
    
    index_split_num = len(xyz_total_copy)/cluster_num
    starting_points = []
    names_of_indices = []
    split_index_range = range(index_split_num, len(xyz_total_copy), index_split_num)
    for index_splitters in split_index_range:
        starting_points.append(xyz_total_copy[index_splitters])
        names_of_indices.append(pdb_obj.get_name_from_xyz(xyz_total_copy[index_splitters]))
    
    #Want to update the starting points by averaging the points in the cluster 5 times
    previous_cluster = []
    while [list(x) for x in previous_cluster] != [list(x) for x in starting_points[:]]:
        #reassign the previous_cluster with current starting points
        previous_cluster = starting_points[:]
        cluster_members = {}
        for xyz in xyz_total_copy:
            distances = [] 
            for xyz2 in starting_points:
                #compare distances and select the smallest distance to assign the xyz1 to that identity
                distance = get_mag_difference(xyz,xyz2)
                distances.append(distance)
            index = distances.index(min(distances))
            try: 
                if cluster_members[index]:
                    cluster_members[index].append(xyz[:])
            except KeyError:
                cluster_members[index] = [xyz[:]]
        #need to average points for a given index and change the starting point
        for cluster_index, xyz_group in cluster_members.iteritems():
            starting_points[cluster_index] = get_geo_center(xyz_group)

    #Here we need to take our points and return a new pdb with the points
    #if chains False, can write to pdb right away otherwise need to calculate differecnes
    if by_chain == False:
        pdb_obj.get_xyz_from_name_list(names_of_indices, update = True)
        pdb_obj.xyz_total = [list(x) for x in starting_points]
    else:
        new_indices = []
        for chain_num in xrange(len(chain_list)):
            new_indices += [x + chain_num*split_index_range[-1] for x in split_index_range]
        pop_counter = 0
        for index in range(len(pdb_obj.xyz_total)):
            if index not in new_indices:
                (pdb_obj.xyz_total).pop(index-pop_counter)
                (pdb_obj.xyz_names).pop(index-pop_counter)
                (pdb_obj.pdb).pop(index-pop_counter)
                pop_counter += 1
            if pdb_obj.xyz_names[index-pop_counter] in names_of_indices:
                pdb_obj.xyz_total[index-pop_counter] = list(starting_points[index-pop_counter])
        
        if pdb_obj.subunits > 1:
            pdb_obj.apply_symmetric_transform_to_main_chain()
    
    for index, xyz in enumerate(pdb_obj.xyz_total):
        pdb_obj.xyz_total[index] = add_vector_difference(xyz[:], orig_center[:])

    out_pdb = pdb_obj.return_to_pdb()
    return out_pdb
    

#### 3.PDB_comparisons ###################################
##
##  Author: Williiam Hansen


#This definition should be used to check the geometry of new loops
#This definition works best determining stretches of atleast 4 residues
#seq_range takes list or string arguements
#It is not necessary but adding a few residues before or after is welcome
def geometric_peptide_angle_check(pdb, seq_range=''):

    pdb = clean_pdb(pdb)

    if type(seq_range) == str:
        seq_range = seq_range.split('-')
    if type(seq_range) == list or type(seq_range) == tuple:
        if len(seq_range) == 2:
            if int(seq_range[0]) < int(seq_range[1]):
                frag = [x for x in pdb if int(x[22:26]) in xrange(int(seq_range[0]), int(seq_range[1]))]
            else:
                frag = [x for x in pdb if int(x[22:26]) in xrange(int(seq_range[1]), int(seq_range[0]))]
        elif len(seq_range) >= 4:
            frag = [x for x in pdb if int(x[22:26]) in seq_range]
        else:
            #print "Did not provide a range of residues or a list of residues greater than 3 residues"
            sys.exit()
    else:
        #print "seq_range must be either str(#-#) or list[#,#]"
        sys.exit()

    #Need to check the information for the given fragment
    #ten pass_checks per residue

    pdb_res_blocks = block_pdb_by_res(frag)

    N_coords = [atom_coordinate(x, 'N', int(x[0][22:26]))[:-1] for x in pdb_res_blocks]
    CA_coords = [atom_coordinate(x, 'CA', int(x[0][22:26]))[:-1] for x in pdb_res_blocks]
    C_coords = [atom_coordinate(x, 'C', int(x[0][22:26]))[:-1] for x in pdb_res_blocks]

    geo_fail_count = 0
    
    def check_angle(atom1, atom2, atom3, floor_val, ceiling_val):
        fail = 0
        try: 
            ang = math.degrees(get_angle(atom1, atom2, atom3))
            if float(ang) < floor_val or float(ang) > ceiling_val:
                fail += 1
        except IndexError:
            return fail
        return fail
        

    for i in xrange(len(pdb_res_blocks)):
        #C-term angle A
        geo_fail_count += check_angle(CA_coords[i], C_coords[i], CA_coords[i+1], 141.5, 150.5)
        geo_fail_count += check_angle(CA_coords[i], C_coords[i], CA_coords[i+2], 115.0, 155.0)
        geo_fail_count += check_angle(CA_coords[i], C_coords[i], CA_coords[i+3], 90.0, 160.0)
        #N-term angle A
        geo_fail_count += check_angle(CA_coords[i], N_coords[i+1], CA_coords[i+1], 152.0, 159.0)
        geo_fail_count += check_angle(CA_coords[i], N_coords[i+2], CA_coords[i+2], 125.0, 150.0)
        geo_fail_count += check_angle(CA_coords[i], N_coords[i+3], CA_coords[i+3], 100.0, 160.0)
        #C-term angle B
        geo_fail_count += check_angle(C_coords[i], CA_coords[i+1], CA_coords[i+2], 84.0, 155.0)
        geo_fail_count += check_angle(C_coords[i], CA_coords[i+2], CA_coords[i+3], 62.0, 160.0)
        #N-term angle B
        geo_fail_count += check_angle(CA_coords[i], CA_coords[i+1], N_coords[i+2], 80.0, 155.0)
        geo_fail_count += check_angle(CA_coords[i], CA_coords[i+1], N_coords[i+3], 60.0, 160.0)
    
    if geo_fail_count == 0:
        return False
    else:
        #print "%s failures were predicted" % geo_fail_count
    #parser.add_argument('-f', '--flip', action="count", default=0, \
        return True

#This is LEGACY, really bad clash check
def clash_check(test_obj, scaffold_obj, mag_check=2.8, exclude_res=[]):
    for i,line in enumerate(test_obj.xyz_total):
        for line2 in scaffold_obj.xyz_total:
            diff = [ line[0] - line2[0], \
                     line[1] - line2[1], \
                     line[2] - line2[2] ]
            mag =  (diff[0]**2)+(diff[1]**2)+(diff[2]**2)
            if mag < mag_check**2:
                return True #if clash se true
            else:
                continue
    #no clashes found
    return False

def gridify_xyz(xyz, grid_size=0.0):
    
    def roundup(x):
        return round(x/math.ceil(grid_size), 0)*math.ceil(grid_size)
    
    low_res = [ roundup(float(xyz[0])),\
                roundup(float(xyz[1])),\
                roundup(float(xyz[2])) ]
    return low_res

#### 4.Vector_math ###################################
##
## Authors: William Hansen, Elliott Dolan, and 
##          Senderle from StackOverflow (Thanks!)

def get_angle(xyz, xyz1, xyz2):
    xyz = numpy.array(xyz) - numpy.array(xyz1)
    xyz2 = numpy.array(xyz2) - numpy.array(xyz1)
    angle = vector_angle(xyz, xyz2)
    return angle
    '''
    xyz = [x - xyz1[i] for i,x in enumerate(xyz)]
    xyz2 = [x - xyz1[i] for i,x in enumerate(xyz2)]
    angle = vector_angle(xyz, xyz2)
    return angle
    '''

def vector_angle2(xyz, xyz2): #LEGACY

    radius_xyz = math.sqrt((xyz[0]**2)+(xyz[1]**2)+(xyz[2]**2))
    radius_xyz2 = math.sqrt((xyz2[0]**2)+(xyz2[1]**2)+(xyz2[2]**2))
    
    dot_over_mag = ((xyz[0]*xyz2[0]) + (xyz[1]*xyz2[1]) + (xyz[2]*xyz2[2])) / (radius_xyz * radius_xyz2)
    if dot_over_mag > 0.9999999 and dot_over_mag < 1.0000001: 
        if get_mag(add_vector_difference(normalize(xyz), normalize(xyz2))) < 0.0001:
            return math.pi
        else:
            return 0.0
            
    elif dot_over_mag < -0.9999999 and dot_over_mag > -1.0000001: 
        if get_mag(add_vector_difference(normalize(xyz), normalize(xyz2))) < 0.0001:
            return math.pi
        else:
            return 0.0
    else:
        theta_angle = math.acos(dot_over_mag)
    return theta_angle


def get_unit_vector(vector):
    return vector / numpy.linalg.norm(vector)

def vector_angle(v1, v2):
    v1_u = get_unit_vector(v1)
    v2_u = get_unit_vector(v2)
    return numpy.arccos(numpy.clip(numpy.dot(v1_u, v2_u), -1.0, 1.0))

def get_mag(xyz):
    x = numpy.array(xyz)
    return numpy.linalg.norm(x)
    '''
    #z = 0.
    #for x in xyz:
        #z += float(x) ** 2.0
    z = sum([float(x)**2 for x in xyz])
    sqZ = math.sqrt(z)
    return sqZ
    '''

def get_mag_difference(xyz1, xyz2):
    xyz = get_vector_difference(xyz1, xyz2)
    return get_mag(xyz)
    '''
    mag = ((xyz1[0] - xyz2[0])**2 +\
           (xyz1[1] - xyz2[1])**2 +\
           (xyz1[2] - xyz2[2])**2)**0.5
   
    return mag
    '''
def get_unit_vector2(xyz): #LEGACY
    xyz_mag = get_mag(xyz)
    if int(xyz_mag) == 0:
        xyz_mag = 1.0
    if type(xyz) == list:
        return [x/xyz_mag for x in xyz]
    #assuming xyz is an numpy array
    else:
        return xyz/xyz_mag

def get_vector_difference(xyz1, xyz2):
    return numpy.array(xyz1) - numpy.array(xyz2)
    '''
    return [xyz1[0]-xyz2[0],\
            xyz1[1]-xyz2[1],\
            xyz1[2]-xyz2[2]]
    '''

def add_vector_difference(xyz1, xyz2):
    return numpy.array(xyz1) + numpy.array(xyz2)
    '''
    return [xyz1[0]+xyz2[0],\
            xyz1[1]+xyz2[1],\
            xyz1[2]+xyz2[2]]
    '''

def get_dihedral(xyz1, xyz2, xyz3, xyz4):
    #if xyz1 == xyz2:
        #print 'dihedral failed'
        #print xyz1
        #print xyz2
        #sys.exit()
    v1 = get_vector_difference(xyz1, xyz2)
    v2 = get_vector_difference(xyz2, xyz3)
    v3 = get_vector_difference(xyz3, xyz4)

    n1 = numpy.cross(v1, v2)
    n2 = numpy.cross(v2, v3)
    
    dot_over_mag = numpy.dot(n1,n2)/(get_mag(n1)*get_mag(n2))
    if round(dot_over_mag, 7) == 1.0:
        return 0.0
    dihedral = math.acos(dot_over_mag)

    #This part is to figure in directionality (- or +)
    m1 = numpy.cross(n1, v2)
    atan_dihedral = math.atan2(numpy.dot(m1,n2), numpy.dot(n1,n2))
    if atan_dihedral < 0:
        dihedral *= -1.0

    return dihedral

#xyz total must be list of xyz lists
def get_geo_center(xyz_total):
    xyz_coord = numpy.array(xyz_total, float)
    xyz_center = numpy.mean(xyz_coord, 0)
    return xyz_center
    
def normalize(xyz, tolerance =0.00001):
    new_xyz = get_unit_vector(xyz)
    return tuple(new_xyz)
    '''
    mag2 = sum(n * n for n in xyz)
    if abs(mag2 - 1.0) > tolerance:
        mag = math.sqrt(mag2)
        if int(mag) == 0:
            mag = 1.0
        xyz = tuple(n / mag for n in xyz)
    return tuple(xyz)
    '''

def cross(xyz, xyz2):
    return numpy.cross(xyz, xyz2)
    '''
    xyz3 = [((xyz[1] * xyz2[2]) - (xyz[2] * xyz2[1])), \
            ((xyz[2] * xyz2[0]) - (xyz[0] * xyz2[2])), \
            ((xyz[0] * xyz2[1]) - (xyz[1] * xyz2[0]))]

    return xyz3
    '''

def dot(xyz, xyz2):
    return numpy.dot(xyz,xyz2)
    '''
    return xyz[0]*xyz2[0] +\
           xyz[1]*xyz2[1] +\
           xyz[2]*xyz2[2]
    '''

def q_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = (w1 * w2) - (x1 * x2) - (y1 * y2) - (z1 * z2)
    x = (w1 * x2) + (x1 * w2) + (y1 * z2) - (z1 * y2)
    y = (w1 * y2) + (y1 * w2) + (z1 * x2) - (x1 * z2)
    z = (w1 * z2) + (z1 * w2) + (x1 * y2) - (y1 * x2)
    return w, x, y, z

def q_conjugate(q):
    q = normalize(q)
    w, x, y, z = q
    return (w, -x, -y, -z)

def q_xyz_mult(q1, xyz):
    v = normalize(xyz)
    p = (0.0, ) + v
    q2 = q_conjugate(q1)
    return q_mult(q_mult(q1, p), q2)[1:]

def axisangle_to_q(v, theta):
    v = normalize(v)
    x, y, z, = v
    theta /= 2
    w = math.cos(theta)
    x = x * math.sin(theta)
    y = y * math.sin(theta)
    z = z * math.sin(theta)
    return w, x, y, z

#### 5.Transform_classes ###################################
##
##  Author: William Hansen

class Transform(object):    

    def __init__(self, pdb):
        if type(pdb) == list or type(pdb) == tuple:    
            ##print type(pose[0])
            self.pdb = pdb
            if type(pdb[0]) == str:
                #xyz_names is resnum_chain_restype_atomtype
                self.xyz_total, \
                self.xyz_names, \
                self.geo_center_mags = pdb_to_xyz(pdb)
                self.grid_name_dict = {}
                self.symm_type = ''
                self.major_symm_axis = ''
                self.subunits = 0
                self.matrix_transform = {}
                self.chain_names = []
                self.sub_symm_objs = {}
                self.sub_symm_centers = {}
        
        #working on this
        '''else:
            self.pose = pose
            ##print type(pose[0])
            self.xyz_total = pose_to_xyz(pose)'''

    ############################################################################
    ##    THe Get Info Definitions
    ## 
    ## Definitions in this set will provide the user with pdb specific
    ## information. The user can update the object in many cases as desired
    ############################################################################
    def pair_name_and_xyz(self):
        return dict(zip(self.xyz_names, self.xyz_total))

    def pair_xyz_and_name(self):
        return dict(zip([str(x) for x in self.xyz_total], self.xyz_names))

    def get_xyz_info(self, num):
        name = self.xyz_names[num]
        coords = self.xyz_total[num]
        #print "Atom index", num, "is named", name, "and is found at", coords
        return [coords, name]

    def get_name_from_xyz(self, xyz):
        xyz = list(xyz)
        for i, val in enumerate(self.xyz_total):
            if xyz == list(val):
                return self.xyz_names[i]
    
    def get_xyz_from_name_list(self, name_list, update=False):
        if type(name_list) == str:
            name_list = [name_list.upper()]
        elif type(name_list) == list:
            name_list = [str(x).upper() for x in name_list]
        coords_by_name = []
        pdb_lines_by_name = []
        mags_by_name = []
        names_by_name = []
        for i, name in enumerate(self.xyz_names):
            if name in name_list:
                coords_by_name.append(self.xyz_total[i])
                pdb_lines_by_name.append(self.pdb[i])
                mags_by_name.append(self.geo_center_mags[i])
                names_by_name.append(name)
        if update == True:
            self.xyz_total = coords_by_name
            self.xyz_names = names_by_name
            self.pdb = pdb_lines_by_name
            self.geo_center_mags = mags_by_name
        return coords_by_name

    def update_coordinates_from_truncated_obj(self, trunc_obj):
        atom_name_list = trunc_obj.xyz_names[:]
        for trunc_name in atom_name_list:
            for i, name in enumerate(self.xyz_names):
                if name == trunc_name:
                    self.xyz_total[i] = trunc_obj.get_xyz_from_name_list(trunc_name)[0][:]
                    mag = get_mag_difference(self.xyz_total[i][:], self.get_geo_center())
                    self.geo_center_mags[i] = mag
    
    #please only merge symmetric things for your sanity and ours
    #with new apply_matrix you can append to main chain and symmeterize the pdb
    #need to code in the append_to_main_chain
    def merge_secondary_object(self, obj2, append_to_main_chain=''):
        for i, name in enumerate(obj2.xyz_names):
            self.xyz_total.append(self.xyz_total[i][:])
            self.xyz_names.append(name)
            self.pdb.append(obj2.pdb[i])
            self.geo_center_mags.append(obj2.geo_center_mags[i])

    def get_geo_center(self):
        xyz_coord = numpy.array(self.xyz_total[:], float)
        xyz_center = numpy.mean(xyz_coord, 0)
        return xyz_center
    
    #if bb set to true 
    def get_symmetric_axes(self, bb=False, normalize_vec=False, vec_size=10.0):
        
        def add_trans(dictionary, translation):
            new_dict = {}
            for k,v in dictionary.iteritems():
                new_v = add_vector_difference(v, translation)
                new_dict[k] = new_v
            return new_dict

        def multiply_mag(dictionary):
            new_dict = {}
            for k,v in dictionary.iteritems():
                v = normalize(v)
                new_v = [vec_size* x for x in v]
                new_dict[k] = new_v
            return new_dict

        def find_major_axis(dictionary):
            dict_keys = dictionary.keys()
            dict_keys.sort()
            dict_keys.sort(key=len, reverse=True)
            self.major_symm_axis = dict_keys[0]
            return 

        chains = set()
        axes = {}
        chain_centers = {}
        main_copy_obj = copy_transform_object(self)
        translation_diff = main_copy_obj.get_geo_center()[:]
        if bb == True:
            main_copy_obj.get_xyz_by_atomtype(['N','CA','C','O'], update=True)
        main_copy_obj.translate([0., 0., 0.])
        oligo_center = main_copy_obj.get_geo_center()
        for name in self.xyz_names:
            chains.add(name.split('_')[1])
        if len(chains) == 1:
            self.symm_type = 'A'
            return [False, axes] #This means that it is a monomer or chainID is same
        #chain_objects = {}
        for chain in chains:
            chain_obj = copy_transform_object(main_copy_obj)
            chain_obj.get_xyz_by_chain(chain, True)
            center = chain_obj.get_geo_center()
            axes[chain] = center
            chain_centers[chain] = center
            #chain_objects[chain] = chain_obj

        self.subunits = len(chains)
        if len(chains) == 2:
            by_chain = zip(*(iter(main_copy_obj.xyz_total[:]),)*(len(main_copy_obj.xyz_total[:])/len(chains)))
            atom_sets_dimer = zip(by_chain[0], by_chain[1])
            define_majors = {}
            for atom_pair in atom_sets_dimer:
                pair_center = numpy.mean(numpy.array(atom_pair), 0)
                
                key = str([int(round(x)) for x in get_unit_vector(list(pair_center))])
                if atom_pair == atom_sets_dimer[0]:
                    define_majors['M1'] = key
                if key not in define_majors.keys():
                    define_majors[key] = [pair_center]
                else:
                    current_dict_items = define_majors[key]
                    current_dict_items.append(pair_center)
                    define_majors[key] = current_dict_items
            for center_name, center_vectors in define_majors.iteritems():
                if center_name == define_majors['M1']:
                    average_vector_center = numpy.mean(numpy.array(center_vectors), 0)
                    chains = list(chains)
                    chains.sort()
                    major_name = '_'.join(chains)
                    chains.sort(reverse=True)
                    minor_name = '_'.join(chains)
                    axes[major_name] = average_vector_center
                    axes[minor_name] = numpy.array(-1. * average_vector_center, float)
            
            if normalize_vec == True:
                axes = multiply_mag(axes)
            axes = add_trans(axes, translation_diff)
            self.symm_type = 'C_2'
            find_major_axis(axes)
            return axes
                
        avg_main_ortho = []
        avg_inverse_ortho = []

        in_plane_monos = {}
        paired_chains = {}
        for key1, val1 in chain_centers.iteritems():
            for key2, val2 in chain_centers.iteritems():
                
                #theta = int(math.degrees(get_angle(val1, oligo_center, val2)))
                theta = round(math.degrees(vector_angle(val1, val2)))
                if theta == 0 or theta == 180:
                    continue
                if theta <= 179 or theta >= 181:
                    real_vector = numpy.cross(val1, val2)
                    ortho_vector = [int(x) for x in real_vector]
                    positives_OV = sum(1 for number in ortho_vector if number > 0)
                    inverse_ortho_vector = [-1*x for x in ortho_vector]
                    positives_IOV = sum(1 for number in inverse_ortho_vector if number > 0)

                    if positives_OV > positives_IOV:
                        main_ortho = ortho_vector
                    else:
                        main_ortho = inverse_ortho_vector
                    
                    unit_real_vector = get_unit_vector(real_vector)
                    possible_majors = [unit_real_vector, -1.*unit_real_vector]
                    
                    if str(main_ortho) in in_plane_monos.keys():
                        key_value = in_plane_monos[str(main_ortho)]
                        if [key2, key1] not in key_value:
                            key_value.append([key1, key2])
                            in_plane_monos[str(main_ortho)] = key_value
                    else:
                        if str([-1*x for x in main_ortho]) not in in_plane_monos.keys():
                            in_plane_monos[str(main_ortho)] = [[key1, key2]]
                    if str(int(round(theta/1.0)*1.0)) in paired_chains.keys():
                        key_theta_value = paired_chains[str(int(round(theta/1.0)*1.0))]
                        if [key2, key1] not in key_theta_value:
                            key_theta_value.append([key1, key2])
                            paired_chains[str(int(round(theta/1.0)*1.0))] = key_theta_value
                    else:
                        paired_chains[str(int(round(theta/1.0)*1.0))] = [[key1, key2]]
        
        if len(in_plane_monos.keys()) == 1: #This asks if all monomers have the same major symmetric axis (Cn)
            #this is a Cn geometry and we need to just grab the major_axis
           
            #make a CA obj and pull xyz, use zip to grab pairs
            broken_into_chains = zip(*(iter(main_copy_obj.xyz_total[:]),)*(len(main_copy_obj.xyz_total[:])/len(chains)))
            atom_sets_by_chain = zip(*broken_into_chains)
            
            for atom_set in atom_sets_by_chain:
                atom_pair_center = numpy.mean(numpy.array(atom_set), 0)
                if list(atom_pair_center) != list(oligo_center):
                    #found first atom_set above or below the plane
                    ang_compare_vals = [ vector_angle(atom_pair_center, possible_majors[0]),\
                                         vector_angle(atom_pair_center, possible_majors[1]) ]
                    smallest_ang_ind = min(xrange(len(ang_compare_vals)),key=ang_compare_vals.__getitem__)
                    break
                
            #symmetric_vector = numpy.array(add_vector_difference(atom_pair_center, translation_diff))
            chains = list(chains)
            chains.sort()
            major_name = '_'.join(chains)
            chains.sort(reverse=True)
            minor_name = '_'.join(chains)
            axes[major_name] = possible_majors[smallest_ang_ind]
            possible_majors.pop(smallest_ang_ind)
            axes[minor_name] = possible_majors[0]
            #axes['minor'] = numpy.array(-1. * atom_pair_center, float)
            
            if normalize_vec == True:
                axes = multiply_mag(axes)
            axes = add_trans(axes, translation_diff)
            self.symm_type = 'C_%s' % len(chains)
            find_major_axis(axes)
            return axes
        else:
            #this is where the Dn axes are determined

            #if the object is 4chains, it is d2 and requires special setup
            #no major or minor defined here
            if len(chains) == 4:
                chain_groups = []
                used_letters = []
                for chain_letter in chains:
                    for chain_letter2 in chains:
                        if chain_letter2 not in used_letters:
                            if chain_letter != chain_letter2:
                                ###you are here doing things with letters
                                chain_groups.append([chain_letter, chain_letter2])
                    used_letters.append(chain_letter)
                for chain_pair in chain_groups:
                    chain_pair.sort()
                    pair_vector = numpy.mean(numpy.array([axes[chain_pair[0]], axes[chain_pair[1]]], float), 0)
                    axes['%s_%s' % (chain_pair[0], chain_pair[1])] = pair_vector
                
                if normalize_vec == True:
                    axes = multiply_mag(axes)            
                axes = add_trans(axes, translation_diff)
                self.symm_type = 'D_2'
                find_major_axis(axes)
                return axes            
                         
            c2_angle = str(min([int(x) for x in paired_chains.keys()]))
            
            largest_list_key = max(paired_chains, key=lambda k: len(paired_chains[k]))
            ##print c2_anglei
            
            top_major = set()
            bottom_major = set()
            for angle, chain_groups in paired_chains.iteritems():
                if angle == largest_list_key:
                    for chain_pair in chain_groups:
                        if chain_pair == chain_groups[0]:
                            top_major.add(chain_pair[0])
                            top_major.add(chain_pair[1])
                        else:
                            if chain_pair[0] in top_major or chain_pair[1] in top_major:
                                top_major.add(chain_pair[0])
                                top_major.add(chain_pair[1])
                            else:
                                bottom_major.add(chain_pair[0])
                                bottom_major.add(chain_pair[1])
                else:
                    #These are your c2 axes store them all
                    for chain_pair in chain_groups:
                        chain_pair.sort()
                        pair_vector = numpy.mean(numpy.array([axes[chain_pair[0]], axes[chain_pair[1]]], float), 0)
                        axes['%s_%s' % (chain_pair[0], chain_pair[1])] = pair_vector
                        
            top_major = list(top_major)
            top_major.sort()
            bottom_major = list(bottom_major)
            bottom_major.sort()
            #symmetric_vector = numpy.mean(numpy.array([axes[x] for x in top_major], float), 0)
            #axes['_'.join(top_major)] = symmetric_vector 
            axes['_'.join(top_major)] = numpy.mean(numpy.array([axes[x] for x in top_major], float), 0)
            #axes['_'.join(bottom_major)] = numpy.array([x * -1. for x in symmetric_vector])
            axes['_'.join(bottom_major)] = numpy.mean(numpy.array([axes[x] for x in bottom_major], float), 0)
            
            #we need to do everything at center and add the translation after

            if normalize_vec == True:
                axes = multiply_mag(axes)            
            axes = add_trans(axes, translation_diff)
            self.symm_type = 'D_%s' % len(top_major)
            find_major_axis(axes)
            return axes            

    def get_main_symmetric_axis(self, d2_main=''):
        copy_obj = copy_transform_object(self)
        axes = copy_obj.get_symmetric_axes()
        split_symm_type = (copy_obj.symm_type).split('_')
        for name in axes.keys():
            if len(name.split('_')) == int(split_symm_type[1]) and 'A' in name.split('_'):
                break
        if split_symm_type == ['D','2'] and d2_main != '':
            return d2_main
        
        return name

    #this can take list or string as arguement
    def get_xyz_by_chain(self, chain, update=False):
        if type(chain) == str:
            chain = [chain.upper()]
        elif type(chain) == list:
            chain = [str(x).upper() for x in chain]
        coords_by_chain = []
        names_by_chain = []
        pdb_lines_by_chain = []
        mags_by_chain = []
        for i, name in enumerate(self.xyz_names):
            if name.split('_')[1] in chain:
                coords_by_chain.append(self.xyz_total[i])
                names_by_chain.append(name)
                pdb_lines_by_chain.append(self.pdb[i])
                mags_by_chain.append(self.geo_center_mags[i])
        if update == True:
            self.xyz_total = coords_by_chain
            self.xyz_names = names_by_chain
            self.pdb = pdb_lines_by_chain
            self.geo_center_mags = mags_by_chain
        return coords_by_chain

    #this can take list or string as arguement
    def get_xyz_by_atomtype(self, atom_type, update=False, reverse=False):
        if type(atom_type) == str:
            atom_type = [atom_type.upper()]
        elif type(atom_type) == list:
            atom_type = [str(x).upper() for x in atom_type]
        coords_by_atom_type = []
        names_by_atom_type = []
        pdb_lines_by_atom_type = []
        mags_by_atom_type = []
        for i, name in enumerate(self.xyz_names):
            if reverse == False:
                if name.split('_')[3] in atom_type:
                    coords_by_atom_type.append(self.xyz_total[i])
                    names_by_atom_type.append(name)
                    pdb_lines_by_atom_type.append(self.pdb[i])
                    mags_by_atom_type.append(self.geo_center_mags[i])
            else:
                if name.split('_')[3] not in atom_type:
                    coords_by_atom_type.append(self.xyz_total[i])
                    names_by_atom_type.append(name)
                    pdb_lines_by_atom_type.append(self.pdb[i])
                    mags_by_atom_type.append(self.geo_center_mags[i])
        if update == True:
            self.xyz_total = coords_by_atom_type
            self.xyz_names = names_by_atom_type
            self.pdb = pdb_lines_by_atom_type
            self.geo_center_mags = mags_by_atom_type
        return coords_by_atom_type
    
    def get_xyz_by_restype(self, res_type, update=False):
        if type(res_type) == str:
            res_type = [res_type.upper()]
        elif type(res_type) == list:
            res_type = [str(x).upper() for x in res_type]
        coords_by_res_type = []
        names_by_res_type = []
        pdb_lines_by_res_type = []
        mags_by_res_type = []
        for i, name in enumerate(self.xyz_names):
            if name.split('_')[2] in res_type:
                coords_by_res_type.append(self.xyz_total[i])
                names_by_res_type.append(name)
                pdb_lines_by_res_type.append(self.pdb[i])
                mags_by_res_type.append(self.geo_center_mags[i])
        if update == True:
            self.xyz_total = coords_by_res_type
            self.xyz_names = names_by_res_type
            self.pdb = pdb_lines_by_res_type
            self.geo_center_mags = mags_by_res_type
        return coords_by_res_type

    def get_xyz_by_resnum(self, res_num, update=False):
        if type(res_num) == str:
            res_num = [res_num.upper()]
        elif type(res_num) == list:
            res_num = [str(x).upper() for x in res_num]
        coords_by_res_num = []
        names_by_res_num = []
        pdb_lines_by_res_num = []
        mags_by_res_num = []
        for i, name in enumerate(self.xyz_names):
            if name.split('_')[0] in res_num:
                coords_by_res_num.append(self.xyz_total[i])
                names_by_res_num.append(name)
                pdb_lines_by_res_num.append(self.pdb[i])
                mags_by_res_num.append(self.geo_center_mags[i])
        if update == True:
            self.xyz_total = coords_by_res_num
            self.xyz_names = names_by_res_num
            self.pdb = pdb_lines_by_res_num
            self.geo_center_mags = mags_by_res_num
        return coords_by_res_num

    def get_xyz_by_grid(self, grid, update=False):
        xyz = {}
        small_name_dict = {}
        for val in grid:
            associated_names = self.grid_name_dict[val]
            for associated_name in associated_names:
                small_name_dict[associated_name] = val

        keys = small_name_dict.keys()
        for i,name in enumerate(self.xyz_names):
            if name in keys:
                xyz[name] = self.xyz_total[i]
        
        if update == False:
            return xyz

        coords_by_grid = []
        names_by_grid = []
        pdb_lines_by_grid = []
        mags_by_grid = []
        for i, name in enumerate(self.xyz_names):
            if name in xyz.keys():
                coords_by_grid.append(self.xyz_total[i])
                names_by_grid.append(name)
                pdb_lines_by_grid.append(self.pdb[i])
                mags_by_grid.append(self.geo_center_mags[i])
        self.xyz_total = coords_by_grid
        self.xyz_names = names_by_grid
        self.pdb = pdb_lines_by_grid
        self.geo_center_mags = mags_by_grid
   
    #This function will change the output pdb lines to reflect the 
    #resnum_chain_resname_atomname found in the name_convert dictionary
    #name_convert dictionary should be a dictionary of keys (original
    #xyz name, with values that correspond to new xyz_name)
    def rename_pdblines_from_xyz_names(self, name_convert):
        #get_xyz_info will give [name, coords] with index input
        final_pdb_lines = []
        final_xyz_names = []
        for ind, pdb_line in enumerate(self.pdb):
            xyz, old_xyz_name = self.get_xyz_info(ind)
            try:
                new_xyz_name = name_convert[old_xyz_name]
                final_xyz_names.append(new_xyz_name)
                split_name = new_xyz_name.split('_')
                line_out = ''.join([pdb_line[:12],\
                                   ('{0:^4}').format(split_name[3]),\
                                   ('{0:^5}').format(split_name[2]),\
                                   split_name[1],\
                                   ('{0:>4}').format(split_name[0]),\
                                   pdb_line[26:]])
                final_pdb_lines.append(line_out)

            except KeyError:
                final_xyz_names.append(old_xyz_name)
                final_pdb_lines.append(pdb_line)
        
        self.pdb = final_pdb_lines[:]
        self.xyz_names = final_xyz_names[:]
        

    def generate_sub_symm_objs(self, symm_type=''):
        #want to produce all possible sub symmetries from the original object
        #Look at symm type and output.
        axes = self.get_symmetric_axes()
        copy_keys = axes.keys()[:]
        copy_keys.sort()
        main_chain = copy_keys[0].split('_')[0]
        for axis in axes.keys():
            ea_chain = axis.split('_')
            if main_chain == axis[0]:
                sub_symm_obj = copy_transform_object(self)
                #Use the get_xyz_by_chain for combinations of letters that include A
                sub_symm_obj.get_xyz_by_chain(ea_chain, update=True)
                original_center = sub_symm_obj.get_geo_center()[:]
                #evaluate symmetry of new object
                new_axes = sub_symm_obj.get_symmetric_axes()
                print new_axes
                name = sub_symm_obj.symm_type + '_' + axis
                if not symm_type:
                    #store in the dictionary
                    self.sub_symm_objs[name] = copy_transform_object(sub_symm_obj)
                    self.sub_symm_centers[name] = original_center[:]
                elif symm_type and symm_type == sub_symm_obj.symm_type:
                    #store in the dictionary
                    self.sub_symm_objs[name] = copy_transform_object(sub_symm_obj)
                    self.sub_symm_centers[name] = original_center[:]
                elif symm_type and symm_type != sub_symm_obj.symm_type:
                    #keep going until we find the desired symm_type
                    continue

    
    def calculate_symmetric_transform(self):
        xyz_total_by_chain = {}
        for index, name in enumerate(self.xyz_names):
            chain = name.split('_')[1]
            if chain not in self.chain_names:
                self.chain_names.append(chain)
            try:
                if xyz_total_by_chain[chain]:
                    xyz_total_by_chain[chain].append(self.xyz_total[index])
            except KeyError:
                xyz_total_by_chain[chain] = [self.xyz_total[index][:]]
        
        xyz_by_chain = {}
        for chain, xyz_total in xyz_total_by_chain.iteritems():
            first_33 = get_geo_center(xyz_total[:len(xyz_total)/3])
            geo_cent = get_geo_center(xyz_total[:])
            final_33 = get_geo_center(xyz_total[len(xyz_total)*2/3:])
            if chain != self.chain_names[0]:
                xyz_by_chain[chain] = [first_33[:], geo_cent[:], final_33[:]]
            else:
                main_chain = [first_33[:], geo_cent[:], final_33[:]]
        
        def calc_row(a1,a2,a3,b1,b2,b3,c1,c2,c3):
            r2 = (b2*a1 - b3*a1 + b3*a2 - b1*a2 + b1*a3 - b2*a3) /\
                 (b1*c3 + b2*c1 + b3*c2 - b1*c2 - b2*c3 - b3*c1)
            r1 = (a1 - a2 - c1*r2 + c2*r2) / (b1 - b2)
            r3 = (a1 - b1*r1 - c1*r2)
            return [r1, r2, r3]

        for chain, xyz_set in xyz_by_chain.iteritems():
            row1 = calc_row( xyz_set[0][0]/main_chain[0][2] , xyz_set[1][0]/main_chain[1][2] ,\
                             xyz_set[2][0]/main_chain[2][2] , main_chain[0][0]/main_chain[0][2] ,\
                             main_chain[1][0]/main_chain[1][2] , main_chain[2][0]/main_chain[2][2] ,\
                             main_chain[0][1]/main_chain[0][2] , main_chain[1][1]/main_chain[1][2] ,\
                             main_chain[2][1]/main_chain[2][2] )
            row2 = calc_row( xyz_set[0][1]/main_chain[0][2] , xyz_set[1][1]/main_chain[1][2] ,\
                             xyz_set[2][1]/main_chain[2][2] , main_chain[0][0]/main_chain[0][2] ,\
                             main_chain[1][0]/main_chain[1][2] , main_chain[2][0]/main_chain[2][2] ,\
                             main_chain[0][1]/main_chain[0][2] , main_chain[1][1]/main_chain[1][2] ,\
                             main_chain[2][1]/main_chain[2][2] )
            row3 = calc_row( xyz_set[0][2]/main_chain[0][2] , xyz_set[1][2]/main_chain[1][2] ,\
                             xyz_set[2][2]/main_chain[2][2] , main_chain[0][0]/main_chain[0][2] ,\
                             main_chain[1][0]/main_chain[1][2] , main_chain[2][0]/main_chain[2][2] ,\
                             main_chain[0][1]/main_chain[0][2] , main_chain[1][1]/main_chain[1][2] ,\
                             main_chain[2][1]/main_chain[2][2] )
            self.matrix_transform[chain] = [row1, row2, row3]
        
        return {k:v for k,v in self.matrix_transform.iteritems()}

    ############################################################################
    ##    The rigid body transform definitions
    ## 
    ## Definitions in this set will change the physical coordinates of the 
    ## input structure self.xyz_total. Generally Rigid body transformations
    ############################################################################

    #rotates object around a given bond designated by two atoms.
    #for object rotations at the origin supply atom1 as [0.0, 0.0, 0.0]
    def rotate_object(self, atom1, atom2, theta, update=True):
        axis = [atom2[0] - atom1[0], \
                atom2[1] - atom1[1], \
                atom2[2] - atom1[2]]
        
        copy_obj = copy_transform_object(self)
        copy_obj.translate([0.0, 0.0, 0.0], atom1[:])
        
        rotate_vector = axisangle_to_q(axis[:], theta)
        
        #rounded_atom2 = [round(x, 6) for x in atom2]
        #rounded_unit_axis = get_unit_vector(rounded_atom2)
        rotated_object = []
        for atom in copy_obj.xyz_total:
            '''
            rounded_atom = [round(x, 6) for x in atom]
            rounded_unit_atom = get_unit_vector(rounded_atom)
            rounded_neg = get_unit_vector([-1*x for x in rounded_atom])
            if rounded_atom == [0.000000, 0.000000, 0.000000]:
                ##print 'xyz   = ' + str(atom)
                Rxyz = atom
            
            elif rounded_unit_atom == rounded_unit_axis:
                Rxyz = atom
            elif rounded_neg == rounded_unit_axis:
                Rxyz = atom
            else:'''
            a_mag = get_mag(atom)
            Rxyz = q_xyz_mult(rotate_vector, atom)
            Rxyz = [(Rxyz[0] * a_mag),\
                    (Rxyz[1] * a_mag),\
                    (Rxyz[2] * a_mag)]

            ##print 'Rxyz  = ' + str(Rxyz)
            rotated_object.append(Rxyz)
        
        if update == True:
            self.xyz_total = rotated_object[:]
            self.translate(atom1, [0.0, 0.0, 0.0])
            return
        else:
            copy_obj.xyz_total = rotated_object[:]
            copy_obj.translate(atom1, [0.0, 0.0, 0.0])
            return copy_obj
        
    
    #axis must be an iterable data member
    def rotate_about_axis(self, axis, theta, update=True):
        
        rotate_vector = axisangle_to_q(axis, theta)
        
        #rounded_axis = [round(x, 6) for x in axis]
        #rounded_unit_axis = get_unit_vector(rounded_axis)
        rotated_object = []
        for atom in self.xyz_total:
            #rounded_atom = [round(x, 6) for x in atom]
            #rounded_unit_atom = get_unit_vector(rounded_axis)
            #rounded_neg = get_unit_vector([-1*x for x in rounded_axis])
            mag = get_mag(atom)
            '''if rounded_atom == [0.000000, 0.000000, 0.000000]:
                Rxyz = atom
            
            elif rounded_unit_atom == rounded_unit_axis:
                Rxyz = atom
            elif rounded_neg == rounded_axis:
                Rxyz = atom
            else:'''
            Rxyz = q_xyz_mult(rotate_vector, atom)
            Rxyz = ((Rxyz[0] * mag),\
                    (Rxyz[1] * mag),\
                    (Rxyz[2] * mag))

            rotated_object.append(Rxyz)
        
        if update == True:
            self.xyz_total = rotated_object[:]
        else:
            new_obj = copy_transform_object(self)
            new_obj.xyz_total = rotated_object[:]
            #new_obj.translate(atom1, [0.0, 0.0, 0.0])
            #self.translate(atom1, [0.0, 0.0, 0.0])
            return new_obj
        

    def translate(self, final, initial = [], update=True):
        ##print "final is " + str(final)
        ##print "intial trans point is " + str(initial)
        if not list(initial):
            ##print "using center, initial is not present"
            input_coord = numpy.array(self.xyz_total, float)
            center = numpy.mean(input_coord, 0)
            translation = [final[0] - center[0], final[1] - center[1], final[2] - center[2]]
        else:
            translation = [final[0] - initial[0], final[1] - initial[1], final[2] - initial[2]]
            ##print "using initial, center not needed"
        ##print ("")
        ##print ("Translation center is: " , translation)
        if update == True:
            self.xyz_total = [ [ x[0] + translation[0],\
                                 x[1] + translation[1],\
                                 x[2] + translation[2] ] for x in self.xyz_total]
            
        else:
            new_obj = copy_transform_object(self)
            new_obj.xyz_total = [ [ x[0] + translation[0],\
                                    x[1] + translation[1],\
                                    x[2] + translation[2] ] for x in new_obj.xyz_total]
            return new_obj

      #axis can be a list or string. String is a key name from get_symmetric_axes
      #this translate just takes a magnitude instead of a point
    def translate_along_axis(self, axis, magnitude=0.0, update=True, invert=False):
        if int(magnitude) == 0:
            return
        obj_symm_axes = self.get_symmetric_axes()
        if type(axis) == str:
            axis = obj_symm_axes[axis]
        if type(axis) == list:
            axis = numpy.array(axis)
        if invert == True:
            axis = [-1.0 * unit for unit in axis]
        #print "axis", axis
        #print "geocenter", self.get_geo_center()[:]
        center = self.get_geo_center()
        axis_to_origin = get_vector_difference(axis, center[:])
        axis_mag = get_mag(axis_to_origin[:])
        unit_vector = [x/axis_mag for x in axis_to_origin]
        #print "axis_mag",axis_mag
        #print "unit_vector", unit_vector
        #print "axis_to_origin", axis_to_origin
        
        
        updated_center = add_vector_difference([x*magnitude for x in unit_vector], center[:])
        if update == True:
            self.translate(updated_center[:], update=True)
            print "after translate geocenter", self.get_geo_center()[:]
        
        else:
            out_obj = self.translate(updated_center[:], update=False)
            return out_obj

    def generate_axes(self):
        xyz_coord = numpy.array(self.xyz_total, float)
        xyz_center = numpy.mean(xyz_coord, 0)
        #center each file with geometric center
        delta = xyz_coord - xyz_center

        #return principal axis matrix
        inertia = numpy.dot(delta.transpose(), delta)
        e_values, e_vectors = numpy.linalg.eig(inertia)

        #The axis we are concerned with at this time is the smallest eigan value (3) but writing in 1 and 2 for later
        for i in xrange(len(e_values)):
            if e_values[i] == max(e_values):
                eval1 = e_values[i]
                axis1 = e_vectors[:,i]
            elif e_values[i] == min(e_values):
                eval3 = e_values[i]
                axis3 = e_vectors[:,i]
            else:
                eval2 = e_values[i]
                axis2 = e_vectors[:,i]

        axis_coords = [axis1, axis2, axis3]

        return axis_coords


    '''def get_terminal_residues(self):
        pdb = self.return_to_pdb()
        CA_xyz = []
        for ind, line in enumerate(pdb):
            if line[11:17].strip() == "CA":
                CA_xyz.append(self.xyz[ind])'''


    #obj1 will align to obj2
    #do not do this for d2
    #only use for c2 and d3+
    def major_symmetric_axis_align(self, obj2, flip=False):
        #Get each objects symmetric axes
        obj1_axes = self.get_symmetric_axes()
        axis1 = self.major_symm_axis
        #print "obj1_axes are:", obj1_axes
        obj2_axes = obj2.get_symmetric_axes()
        axis2 = obj2.major_symm_axis
        #print "obj2_axes are:", obj2_axes
        #print obj1_axes
        #print obj2_axes

        if flip == False:
            print axis2, 'will be aligned with', axis1
            if list(get_unit_vector(obj1_axes[axis1])) == list(get_unit_vector(obj2_axes[axis2])):
                return
            ortho_vector = numpy.cross(obj1_axes[axis1][:], obj2_axes[axis2][:])
            theta = vector_angle(obj1_axes[axis1][:], obj2_axes[axis2][:])
            #print 'o_vector', ortho_vector
            #print 'theta', theta
        else:
            print axis2, 'will be aligned with the opposite of', axis1
            flipped = [-1.*x for x in obj1_axes[axis1]]
            if list(get_unit_vector(flipped[:])) == list(get_unit_vector(obj2_axes[axis2])):
                return
            ortho_vector = numpy.cross(flipped[:], obj2_axes[axis2][:])
            theta = vector_angle(flipped[:], obj2_axes[axis2][:])
        
        #this is totally wrong, need to fix 
        if theta > math.pi - 0.1 or theta < 0.1:
            #pre_final = self.get_geo_center()
            #self.rotate_object([0.,0.,0.], 
            #self.translate([0.,10.,0.])
            #obj1_axes_fix = self.get_symmetric_axes()
            self.rotate_about_axis([10.0,0.0,0.0], math.pi/16.0)
            self.rotate_about_axis([0.0,10.0,0.0], math.pi/16.0)
            #self.translate(pre_final[:])
            obj1_axes_new = self.get_symmetric_axes()
            if flip == False:
                ortho_vector = numpy.cross(obj1_axes_new[axis1][:], obj2_axes[axis2][:])
                theta = vector_angle(obj1_axes_new[axis1][:], obj2_axes[axis2][:])
            else:
                flipped = [-1.*x for x in obj1_axes_new[axis1]]
                ortho_vector = numpy.cross(flipped[:], obj2_axes[axis2][:])
                theta = vector_angle(flipped[:], obj2_axes[axis2][:])
        
        #print 'before rotate'
        self.rotate_object([0.0,0.0,0.0], normalize(ortho_vector), theta)
    
    def axis_align(self, obj2, pre_axis_n, post_axis_n, rotate_axis=[]):
        
        obj1_axes = self.get_symmetric_axes()
        obj2_axes = obj2.get_symmetric_axes()
        pre_axis = obj1_axes[pre_axis_n][:]
        post_axis = obj2_axes[post_axis_n][:]

        if list(get_unit_vector(pre_axis)) == list(get_unit_vector(post_axis)):
            return
        ortho_vector = numpy.cross(pre_axis[:], post_axis[:])
        theta = vector_angle(pre_axis[:], post_axis[:])
        
        
        if theta > math.pi - 0.1 or theta < 0.1:
            pre_final = self.get_geo_center()
            self.translate([0.,10.,0.])
            obj1_axes_fix = self.get_symmetric_axes()
            self.rotate_about_axis([10.0,0.0,0.0], math.pi/16.0)
            self.translate(pre_final)
            obj1_axes_new = self.get_symmetric_axes()
            ortho_vector = numpy.cross(obj1_axes_new[pre_axis_n][:], post_axis[:])
            theta = vector_angle(obj1_axes_new[pre_axis_n][:], post_axis[:])
        
        #print 'before rotate'
        if rotate_axis:
            cen_spot = self.get_geo_center()
            self.rotate_object(cen_spot[:], normalize(rotate_axis), theta)
        else:
            cen_spot = self.get_geo_center()
            self.rotate_object(cen_spot[:], normalize(ortho_vector), theta)
            

    #moving axis is axis1, static axis is axis2 
    def align_object_to_vector(self, axis1, axis2, update=True):
        #determin protein once
        if list(get_unit_vector(axis1)) == list(get_unit_vector(axis2)):
            #print "Axis1 and Axis2 in align_object_to_vector cannot be the same value"
            return
        ortho_vector = numpy.cross(axis1, axis2)
        theta = vector_angle(axis1, axis2)
        if theta > math.pi - 0.0001 or theta < 0.0001:
            theta_2 = vector_angle(axis1, [0.0,1.0,0.0])
            ortho_vector = numpy.cross(axis1, [0.0,1.0,0.0])
            if theta_2 > 179.99 or theta_2 < 0.01:
                ortho_vector = numpy.cross(axis1, [1.0,1.0,0.0])

        if update == True:
            self.rotate_object([0.0,0.0,0.0], ortho_vector, theta)
        else:
            new_obj = copy_transform_object(self)
            new_obj.rotate_object([0.0,0.0,0.0], ortho_vector, theta)
            return new_obj

    def get_transform_operations(self):
        copy_obj = copy_transform_object(self)
        copy_obj.get_xyz_by_atomtype(['N','CA','C','O'], True)
        trans = copy_obj.get_geo_center()[:]
        copy_obj.translate([0.,0.,0.])
        
        #copy_axes = copy_obj.get_symmetric_axes()
        #major_axis = copy_axes['major']
        #take first 10 residues and find the geocenter later to make this more robust
        major_axis = copy_obj.xyz_total[0]
        major_rot = [numpy.cross([0.,0.,1.], major_axis), vector_angle([0.,0.,1.], major_axis)]
        copy_obj.align_object_to_vector(major_axis, [0.,0.,1.], True)

        #copy_axes = copy_obj.get_symmetric_axes()
        #A_axis = copy_axes['A']
        A_axis = copy_obj.xyz_total[-1]
        A_rot = [numpy.cross([1.,0.,0.], A_axis), vector_angle([1.,0.,0.], A_axis)]
        
        operations = [A_rot, major_rot, trans]
        return operations
 
    def apply_transform_operations(self, operations, rvs=False):
        self.translate([0.,0.,0.])
        #symm_axes = self.get_symmetric_axes()
        #major_axis = copy_axes['major']
        major_axis = self.xyz_total[0]
        self.align_object_to_vector(major_axis, [0.,0.,1.], True)

        #copy_axes = copy_obj.get_symmetric_axes()
        #A_axis = copy_axes['A']
        A_axis = self.xyz_total[-1]
        self.align_object_to_vector(A_axis, [1.,0.,0.], True)

        if rvs == False:
            #once object is at starting position, we apply operations in rvs order
            self.rotate_object([0.,0.,0.],operations[0][0], operations[0][1], True)
            self.rotate_object([0.,0.,0.],operations[1][0], operations[1][1], True)
            self.translate(operations[2])
        
        else:
            #once object is at starting position, we apply operations in rvs order
            self.rotate_object([0.,0.,0.],operations[0][0], -1.0*operations[0][1], True)
            self.rotate_object([0.,0.,0.],operations[1][0], -1.0*operations[1][1], True)
            self.translate([-1.0*x for x in operations[2]])
    
    def apply_symmetric_transform_to_main_chain(self):
        if not self.matrix_transform:
            print "Cannot apply symmetric transformation, object has no transform"
            sys.exit()
        self.get_xyz_by_chain(self.chain_names[0], update=True)
        main_chain_xyz = self.xyz_total[:]
        for chain in self.chain_names[1:]:
            for index, xyz in enumerate(main_chain_xyz):
                x = (self.matrix_transform[chain][0][0] * xyz[0]) +\
                    (self.matrix_transform[chain][0][1] * xyz[1]) +\
                    (self.matrix_transform[chain][0][2] * xyz[2])
                y = (self.matrix_transform[chain][1][0] * xyz[0]) +\
                    (self.matrix_transform[chain][1][1] * xyz[1]) +\
                    (self.matrix_transform[chain][1][2] * xyz[2])
                z = (self.matrix_transform[chain][2][0] * xyz[0]) +\
                    (self.matrix_transform[chain][2][1] * xyz[1]) +\
                    (self.matrix_transform[chain][2][2] * xyz[2])
                
                self.xyz_total.append([x,y,z])
                equiv_name = self.xyz_names[index].split('_')
                equiv_name[1] = chain
                self.xyz_names.append('_'.join(equiv_name))
                self.pdb.append(self.pdb[index][:21] + chain + self.pdb[index][22:])
        return

    #This definition is different from the above in that it applys a pregenerated symmetric transform
    #pregenerated symmetric transform to the transform object and returns each of the new objects
    def apply_transform_to_object(self, matrix):
        if not matrix:
            print "Cannot apply symmetric transformation, no matrix supplied"
            sys.exit()
        main_chain_xyz = self.xyz_total[:]
        transformed_objects = []
        for chain, xyz_set in matrix.iteritems():
            new_object = copy_transform_object(self)
            for index, xyz in enumerate(main_chain_xyz):
                x = (xyz_set[0][0] * xyz[0]) +\
                    (xyz_set[0][1] * xyz[1]) +\
                    (xyz_set[0][2] * xyz[2])
                y = (xyz_set[1][0] * xyz[0]) +\
                    (xyz_set[1][1] * xyz[1]) +\
                    (xyz_set[1][2] * xyz[2])
                z = (xyz_set[2][0] * xyz[0]) +\
                    (xyz_set[2][1] * xyz[1]) +\
                    (xyz_set[2][2] * xyz[2])
                
                new_object.xyz_total[index] = [x,y,z]
            transformed_objects.append(copy_transform_object(new_object))
        return transformed_objects
    ############################################################################
    ##    The comparison and self-defining definitions
    ## 
    ## Definitions in this set will compare one object to another. Also defs
    ## that learn and build up the object information are found here.
    ############################################################################
    def get_atom_grid(self, grid_size=0.0, atom_name='', refinement=False):
        
        def roundup(x):
            return round(x/math.ceil(grid_size), 0)*math.ceil(grid_size)
            #return int(math.ceil(int(x)/int(math.ceil(grid_size))))*int(math.ceil(grid_size))
        def add_grid_member(grid_name, val):
            if grid_name in self.grid_name_dict.keys():
                current_grid_names = self.grid_name_dict[grid_name]
                current_grid_names.append(val)
            else:
                self.grid_name_dict[grid_name] = [val]
        
        grid = set()
        self.grid_name_dict = {}
        if atom_name:
            xyz_set = self.get_xyz_by_atomtype(atom_name)
        else:
            xyz_set = self.xyz_total[:]
        for x in xyz_set:
            pre_low_res = [ roundup(float(x[0])),\
                            roundup(float(x[1])),\
                            roundup(float(x[2])) ]
            #self.grid_name_dict[self.xyz_names[i]] = low_res
            name = self.get_name_from_xyz(x)
            low_res = ( str(pre_low_res) )
            add_grid_member(low_res, name)
            grid.add(low_res)
            
            if refinement == True:
                extra_low_res = []
                for ind, val in enumerate(pre_low_res):
                    val_copy_plus = pre_low_res[:]
                    val_copy_minus = pre_low_res[:]
                    val_copy_plus[ind] += math.ceil(grid_size)
                    val_copy_minus[ind] += math.ceil(grid_size)
                    low_res_plus = ( str(val_copy_plus) )
                    low_res_minus = ( str(val_copy_minus) )
                    add_grid_member(low_res_plus, name)
                    add_grid_member(low_res_minus, name)
                    grid.add(low_res_plus)
                    grid.add(low_res_minus)

        return grid

    def clash_check(self, test_scaff, mag_check=2.8, \
                    include_atom_type=[], include_resnum=[], \
                    include_chain=[], refine=False):
        
        proto = copy_transform_object(self)
        if include_chain:
            proto.get_xyz_by_chain(include_chain, True)
        if include_resnum:
            proto.get_xyz_by_resnum(include_resnum, True)
        if include_atom_type:
            proto.get_xyz_by_atomtype(include_atom_type, True)
        
        #proto_grid, proto_name_grid = make_clash_grid(proto)
        proto_grid = proto.get_atom_grid(mag_check+1.0, refinement=refine)

        mega = copy_transform_object(test_scaff)
        #mega_grid, mega_name_grid = make_clash_grid(mega)
        mega_grid = mega.get_atom_grid(mag_check+1.0, refinement=refine)
        
        clashing_grids = proto_grid & mega_grid
        #print '\n clashing grids', clashing_grids
        if not clashing_grids:
            #print '\n no grid \n'
            return False
        

        #proto_coords_to_check = return_xyz_coords_to_check(proto, clashing_grids, proto_name_grid)
        proto.get_xyz_by_grid(clashing_grids, True) 
        #mega_coords_to_check = return_xyz_coords_to_check(mega, clashing_grids, mega_name_grid)
        mega.get_xyz_by_grid(clashing_grids, True) 
        
        #print proto.xyz_names
        #print '\n'
        #print mega.xyz_names
        for xyz1 in proto.xyz_total:
            for xyz2 in mega.xyz_total:
                diff = [ xyz2[0] - xyz1[0], \
                         xyz2[1] - xyz1[1], \
                         xyz2[2] - xyz1[2] ]
                mag =  (diff[0]**2)+(diff[1]**2)+(diff[2]**2)
                if mag < mag_check**2:
                    ##print "Atoms", name1, "and", name2, "are clashing.\nDistance =", mag**0.5
                    return True #if clash se true
                else:
                    continue
        #no clashes found
        return False

    def rmsd_check(self, q_obj):
        #fix this later to make it more robust
        xyz_1 = self.xyz_total[:]
        xyz_2 = q_obj.xyz_total[:]
        if len(xyz_1) == len(xyz_2):
            mag_squared_sum = 0
            for i, xyz in enumerate(xyz_1):
                mag_squared_sum += ((xyz[0] - xyz_2[i][0])**2)+\
                                    ((xyz[0] - xyz_2[i][0])**2)+\
                                    ((xyz[0] - xyz_2[i][0])**2)
            rmsd = (mag_squared_sum/len(xyz))**0.5
            return rmsd
        else:
            print 'rmsd not calculated, total atoms not the same'
            return 'NA'

    ############################################################################
    ##    The output definitions
    ## 
    ## All output methods are found in this block. Return an object to pdb 
    ## format or directly write/append an object to file as single or model
    ############################################################################

    def return_to_pdb(self):
        out_pdb = []
        for i, coord in enumerate(self.xyz_total):
            if abs(coord[0]) >= 100. and abs(coord[0]) < 1000. and coord[0] < 0.:
                x_comp = "%.2f" % coord[0]
            elif abs(coord[0]) >= 1000. and abs(coord[0]) < 10000.:
                x_comp = "%.2f" % coord[0]
                if coord[0] < 0.:
                    x_comp = "%.1f" % coord[0]
            elif abs(coord[0]) >= 10000. and abs(coord[0]) < 100000.:
                x_comp = "%.1f" % coord[0]
                if coord[0] < 0.:
                    x_comp = "%.0f" % coord[0]
            elif abs(coord[0]) >= 100000. and abs(coord[0]) < 1000000.:
                x_comp = int(coord[0])
            else:
                x_comp = "%.3f" % coord[0]

            if abs(coord[1]) >= 100. and abs(coord[1]) < 1000. and coord[1] < 0.:
                y_comp = "%.2f" % coord[1]
            elif abs(coord[1]) >= 1000. and abs(coord[1]) < 10000.:
                y_comp = "%.2f" % coord[1]
                if coord[1] < 0.:
                    y_comp = "%.1f" % coord[1]
            elif abs(coord[1]) >= 10000. and abs(coord[1]) < 100000.:
                y_comp = "%.1f" % coord[1]
                if coord[1] < 0.:
                    y_comp = "%.0f" % coord[1]
            elif abs(coord[1]) >= 100000. and abs(coord[1]) < 1000000.:
                y_comp = int(coord[1])
            else:
                y_comp = "%.3f" % coord[1]

            if abs(coord[2]) >= 100. and abs(coord[2]) < 1000. and coord[2] < 0.:
                z_comp = "%.2f" % coord[2]
            elif abs(coord[2]) >= 1000. and abs(coord[2]) < 10000.:
                z_comp = "%.2f" % coord[2]
                if coord[2] < 0.:
                    z_comp = "%.1f" % coord[2]
            elif abs(coord[2]) >= 10000. and abs(coord[2]) < 100000.:
                z_comp = "%.1f" % coord[2]
                if coord[2] < 0.:
                    z_comp = "%.0f" % coord[2]
            elif abs(coord[2]) >= 100000. and abs(coord[2]) < 1000000.:
                z_comp = int(coord[2])
            else:
                z_comp = "%.3f" % coord[2]


            line_out = ''.join([self.pdb[i][0:30],\
                               ('{0:>8}').format(x_comp),\
                               ('{0:>8}').format(y_comp),\
                               ('{0:>8}').format(z_comp),\
                               self.pdb[i][54:]])
            out_pdb.append(line_out)
        return out_pdb

    def write_to_file(self, path, chain=''):
        if not chain:
            with open(path, 'w') as myfile:
                myfile.writelines(self.return_to_pdb())
        else:
            pdb = self.return_to_pdb()
            pdb = [x for x in pdb if str(x[21]) == chain]
            with open(path, 'w') as myfile:
                myfile.writelines(pdb)
            
                
    
    #if any('MODEL' in s for s in old_pdb): might still be used
    #in place of the try except
    def append_to_file(self, path, model=False):
        if model == False:
            with open(path, 'a') as myfile:
                myfile.writelines(self.return_to_pdb())
        else:    
            try:
                old_pdb = read_file(path)

                try:
                    indices = [i for i, s in enumerate(old_pdb) if 'MODEL' in s] 
                    last_model_index = indices[-1]
                    last_model = int(old_pdb[last_model_index].split()[1])
                    #print "%s contains MODELS, appending %s to file." % (path, last_model + 1)
                    output = ['MODEL  %s\n' % (last_model +1)] + self.return_to_pdb() + ['ENDMDL\n']
                except ValueError:
                    #print "Value Error Raised, path submitted does not contain MODELS"
                    output = ['MODEL  1\n'] + old_pdb + ['ENDMDL\n'] + \
                             ['MODEL  2\n'] + self.return_to_pdb() + ['ENDMDL\n']

                    with open(path, 'w') as myfile:
                        myfile.writelines(output)
                    return

            except IOError:
                #print "%s does not exist. Writing as new file." % path
                output = ['MODEL  1\n'] + self.return_to_pdb() + ['ENDMDL\n']
                with open(path, 'w') as myfile:
                    myfile.writelines(output)
                return
                
            with open(path, 'a') as myfile:
                myfile.writelines(output)


class TransformCopy(Transform):    
    
    def __init__(self, info):
        self.pdb = info['pdb']
        self.xyz_total = info['xyz_total']
        self.xyz_names = info['xyz_names']
        self.geo_center_mags = info['geo_center_mags']
        self.grid_name_dict = info['grid_name_dict']
        self.symm_type = info['symm_type']
        self.major_symm_axis = info['major_symm_axis']
        self.subunits = info['subunits']
        self.matrix_transform = info['matrix_transform']
        self.chain_names = info['chain_names']
        self.sub_symm_objs = info['sub_symm_objs']
        self.sub_symm_centers = info['sub_symm_centers']

        super(TransformCopy, self).__init__('')
    
    def assign(self, obj):
        self.__dict__ = obj.__dict__.copy()
        return
        

class TransformVector(Transform):
    
    def __init__(self, xyz):
        self.xyz = xyz
        super(TransformVector, self).__init__(self.xyz)

    #this is assuming axis is coming from origin and must be given args in a for loop
    #for larger objects it is best to use rotate_object and provide atom1 as [0.0,0.0,0.0]
    def rotate_vector(self, axis, theta):

        mag =  ((self.xyz[0]**2)+(self.xyz[1]**2)+(self.xyz[2]**2))**0.5

        rotate_vector = axisangle_to_q(axis, theta)

        Rxyz = q_xyz_mult(rotate_vector, self.xyz)
        self.xyz = ((Rxyz[0] * mag), (Rxyz[1] * mag), (Rxyz[2] * mag))

    def translate(self, final, update=True):
        translation = [final[0] - self.xyz[0], final[1] - self.xyz[1], final[2] - self.xyz[2]]
        if update == True:
            self.xyz  =  [ self.xyz[0] + translation[0],\
                           self.xyz[1] + translation[1],\
                           self.xyz[2] + translation[2] ]
            
        else:
            return   [ self.xyz[0] + translation[0],\
                       self.xyz[1] + translation[1],\
                       self.xyz[2] + translation[2] ]

      #axis can be a list or string. String is a key name from get_symmetric_axes
      #this translate just takes a magnitude instead of a point
    def translate_along_axis(self, axis, magnitude=0.0, update=True, invert=False):
        if int(magnitude) == 0:
            return
        
        if invert == True:
            axis = [ -1.0 * unit for unit in axis]
        
        axis_to_origin = get_vector_difference(axis, self.xyz)
        axis_mag = get_mag(axis_to_origin)
        unit_vector = [x/axis_mag for x in axis_to_origin]
        
        if update == True:
            self.translate(add_vector_difference([x*magnitude for x in unit_vector],\
                       self.xyz[:]), True)
        
        else:
            return self.translate(add_vector_difference([x*magnitude for x in unit_vector],\
                       self.xyz[:]), False)

def read_chi_dist_file(file_path):
    with open(file_path, 'r') as f:
        my_file = f.readlines()
    blocks = [i.split() for i in my_file]
    final_blocks = []
    block = []
    for item in blocks:
        if item == blocks[-1]:
            if not block:
                continue
            else:
                block.append(item)
                final_blocks.append(block)
                block = []
        elif not item:
            if not block:
                continue
            else:
                final_blocks.append(block)
                block = []
            continue
        else:
            block.append(item)
    return final_blocks

def read_SyPRIS_file(file_path):
    with open(file_path, 'r') as f:
        my_file = f.readlines()
        my_file = [x.strip('\r\n').split(',') for x in my_file]
    pdbs = [x[0] for x in my_file]
    residues = [x[1] for x in my_file]
    rigid_rots = [x[2] for x in my_file]
    rmsds = [x[3] for x in my_file]
    anglelogs = [x[4] for x in my_file]
    chis = [x[5:] for x in my_file]
    
    return pdbs, residues, rigid_rots, rmsds, anglelogs, chis

'''
#for angle_log to work it requires a list of vectors
#example: [ [[1,2,3],[4,5,6]], [[],[]] ]
#where paired atom coords would produce a vector
#from atom1 to atom2 threshold default is 20 degrees, 
#and must be supplied as radians
def angleLog(test, scaffold, threshold=(math.pi/9.0)):
    compare_vectors = zip(test, scaffold)
    vector_ang_sum = 0.00001
    for num_vectors, vector_pairs in enumerate(compare_vectors):
        vector1 = Transform(vector_pairs[0])
        vector1.translate([0.0,0.0,0.0], vector_pairs[0][0])
        vector2 = Transform(vector_pairs[1])
        vector2.translate([0.0,0.0,0.0], vector_pairs[1][0])
        vector_ang_sum += vector_angle(vector1.xyz_total[1], vector2.xyz_total[1])
    angle_log = math.log10((vector_ang_sum) / (len(compare_vectors) * threshold))
    return angle_log
'''
