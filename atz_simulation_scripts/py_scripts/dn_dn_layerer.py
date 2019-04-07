#switched a pop statement with a del statement
#added a dictionary that keeps track of which oligomers are associated with the others
#added two values that need to be deleted. The have 'delete' in the comment

#!/usr/bin/env python
from transform import *
from AssemB_pkgs import *
import os, sys, math, numpy, time, csv, optparse, collections

#These values are necessary for the D3 and are hardcoded
incompatible_axes = { 'A_E' : ['A_F', 'B_E', 'A_E'],\
                      'B_D' : ['B_E', 'C_D', 'B_D'],\
                      'C_F' : ['C_D', 'A_F', 'C_F'],\
                      'A_F' : ['A_E', 'C_F', 'A_F'],\
                      'B_E' : ['A_E', 'B_D', 'B_E'],\
                      'C_D' : ['B_D', 'C_F', 'C_D'] }

compatible_axes   = { 'A_E' : ['B_D', 'C_F', 'C_D'],\
                      'B_D' : ['A_E', 'C_F', 'A_F'],\
                      'C_F' : ['A_E', 'B_D', 'B_E'],\
                      'A_F' : ['B_E', 'C_D', 'B_D'],\
                      'B_E' : ['A_F', 'C_D', 'C_F'],\
                      'C_D' : ['A_F', 'B_E', 'A_E'] }

equal_axes        = { 'A_E' : ['B_D', 'C_F'],\
                      'B_D' : ['A_E', 'C_F'],\
                      'C_F' : ['A_E', 'B_D'],\
                      'A_F' : ['B_E', 'C_D'],\
                      'B_E' : ['A_F', 'C_D'],\
                      'C_D' : ['A_F', 'B_E'] }


def swap(var1, var2):
    var_ref = var1
    var1 = var2
    var2 = var_ref
    return var1, var2

def dn_dn_layerer( static_filename, mobile_filename, layers,\
                   trans_bin, trans_prob, rot_bin, rot_prob,\
                   null1, null2, paired_bins, paired_prob,\
                   d2_main_axis, clash_mag, a_frac_sele, b_frac_sele,\
                   file_name, static_axis_bin, mobile_axis_bin ):
    #arb unit, delete
    particle_goal = 5000
    
    init_null = null2
    init_frac = b_frac_sele
    temp_file_name = file_name + '_null' + str(init_null) + '_frac' + str(init_frac) 
    two_comp_obj = Dn_D2_AssemMaker(static_filename, mobile_filename)
    #Each layer needs to keep track of all possible static objects
    #We can keep to reference mobile objects and assign them each time
    #We must alternate between scaffold types while iterating
    two_comp_obj.static_obj.translate([0.,0.,0.],update=True)
    final_scaffold_grid = two_comp_obj.static_obj.get_atom_grid(grid_size=clash_mag)
    list_of_accepted_additions = [copy_transform_object(two_comp_obj.static_obj)]
    
    a_not_last_layer = 1
    a_count = 1
    b_count = 0
    branch_sum = 0.0
    
    distances = []
    rotations = []
    number_neighbors = {}
    
    #arb unit, delete
    layers = 10000
    for layer in xrange(1,layers+1):
        
        if layer == 1:
            #current static_objs are the objs in the layer that we will place a new obj on
            current_static_objs = [ [copy_transform_object(two_comp_obj.static_obj), ''] ]
        else:
            #swap which is mobile and which in static
            two_comp_obj.swap_object_properties()
            #swap nulls, main axes, and fractions
            null1, null2 = swap(null1, null2)
            #static_main, mobile_main = swap(static_main, mobile_main)
            a_frac_sele, b_frac_sele = swap(a_frac_sele, b_frac_sele)

            if paired_bins:
                static_axis_bin, mobile_axis_bin = swap(static_axis_bin, mobile_axis_bin)
            #current_static_objs = {k:v for k,v in last_layer_placed.iteritems()}
            #last layer placed should be all of the accepted pdbs from the previous layer
            current_static_objs = last_layer_placed[:]
            if not current_static_objs:
                #if previous layer had nothing in it we need to break 
                break
        
        mobile_main = two_comp_obj.mobile_obj.get_main_symmetric_axis(d2_main_axis)
        
        #we will set the active frac equal to the fraction of the mobile we want to place this round
        active_frac = b_frac_sele
        print 'layer:', layer

        last_layer_placed = []
        associated_statics = []
        rigid_body_transformation = []
        random.shuffle(current_static_objs)
        #depending on the ratio of oligomers, we may need to reduce the amount
        cut_index = int(len(current_static_objs) * active_frac + 1)
        #1,3,5,7 are the layers when A is static.
        if layer % 2 != 0 and layer != 1: 
            #append to the branch_sum the number of statics in the layer if it is odd
            #this will account for all Cs connected to the A before placement.
            branch_sum += len(current_static_objs)
        #use the cut index to cut the number of statics used to place mobiles
        current_static_objs = current_static_objs[:cut_index]
        #static obj is an obj, static_axis_used is either '' if layer 1 or an axis
        #we know that we can only allow certain axes based on the static_axis_used
        for static_object, static_axis_used in current_static_objs:
             
            static_object.get_symmetric_axes()
            
            stat_geo = static_object.get_geo_center()
            stat_geo_name = str(static_object.subunits) + ',' + str([int(x) for x in stat_geo])
            if stat_geo_name not in number_neighbors.keys():
                number_neighbors[stat_geo_name] = set()
            
            static_main = static_object.get_main_symmetric_axis(d2_main_axis)
            print static_object.get_symmetric_axes()
            print 'static_main', static_main
            print 'mobile_main', mobile_main
            #append_max should be based on static_object number of chains
            if layer == 1:
                #max is three if 1st static
                append_max = (static_object.subunits/2)
            else:
                #max is 1 less than half the chains, 2 for 6 and 1 for 4
                append_max = (static_object.subunits/2) - 1
            #need to pick from the bins randomly.
            #need to determine how to terminate
            #can hardcode for this problem
            termination = False
            append_count = 0
            append_obj = 0
            #these are the static axes used while placing mobiles around one static
            static_axes_used = []
            #want to keep doing this until we have either nulled out or appended max
            key_error_count = 0
            while append_max > 0:
                #make an allowed axes bin
                if static_axis_used == '':
                    
                    axes_not_allowed = []
                else:
                    axes_not_allowed = incompatible_axes[static_axis_used][:] + [static_axis_used]
                if static_axes_used:
                    for axis in static_axes_used:
                        axes_not_allowed += incompatible_axes[axis][:]
                        axes_not_allowed += [axis]


                #check to see if the bin probabilities are paired up
                if paired_bins == True:
                    #zip together the paired bins
                    axis_trans_rot_prob = zip(static_axis_bin, mobile_axis_bin, trans_bin, rot_bin, paired_prob)
                    rand_frac = random.random()
                    #zip together probs with some index that we use to access the axis_trans_rot_prob
                    tuple_probs = zip(paired_prob, range(len(paired_prob)))
                    #should sort by the probs
                    tuple_probs.sort()
                    tuple_probs.reverse()
                    sum_of_probs = 0.0
                    for prob, ind in tuple_probs:
                        #go through each prob (most likely to least) if your rand frac is less than your prob, you use it
                        if prob + sum_of_probs > rand_frac:
                            paired_values = axis_trans_rot_prob[ind][:]
                            #possible static_axes
                            #HARDCODED
                            if static_main == 'A_D':
                                static_axes = ['A_C']
                                possible_mobile_axes = equal_axes[paired_values[1]][:] + [paired_values[1]]
                                random.shuffle(possible_mobile_axes)
                                mobile_axes = [possible_mobile_axes[0]]
                                rotation = 360.0-paired_values[3]
                            else:
                                #static axes will be the static axis and all axes similar to the one in the axis bin
                                possible_static_axes = equal_axes[paired_values[0]][:] + [paired_values[0]]
                                #remove static axes that we have used already
                                static_axes = [axis for axis in possible_static_axes if axis not in axes_not_allowed]
                                if not static_axes:
                                    sum_of_probs += prob
                                    continue
                                mobile_axes = [paired_values[1]]
                                rotation = paired_values[3]
                            #shuffle the static axes and choose randomly
                            random.shuffle(static_axes)
                            translation = paired_values[2]
                            break
                        sum_of_probs += prob
                

                if not static_axes:
                    key_error_count += 1
                    if key_error_count < 3:
                        continue
                    else:
                        break
                #static_axes are passed to here from the previous conditionals
                two_comp_obj.static_obj = copy_transform_object(static_object)
                print static_axes[0]
                print mobile_axes[0]
                print static_main
                print mobile_main
                print translation
                print rotation
                new_mobile_obj = two_comp_obj.add_mobile_to_static( translation, rotation, null2,\
                                                                      static_axes[0], mobile_axes[0],\
                                                                      static_main, mobile_main )
                if new_mobile_obj == None:
                    print 'no objs in layer key error, appendmax stays the same'
                    key_error_count += 1
                    if key_error_count < 3:
                        continue
                    else:
                        break
                if not new_mobile_obj:
                    append_max -=1
                    print 'no objs in layer hit null, append count is:', append_max
                    continue
                last_layer_placed += [[new_mobile_obj, mobile_axes[0]]]
                static_axes_used.append(static_axes[0])
                rigid_body_transformation.append((translation, rotation))

                #mobile_geo name is used to keep track of each oligomer in a dictionary
                #It will store the type of oligomer and the geo ceneter it is located at as well as # of neighbors
                mobile_geo = new_mobile_obj.get_geo_center()
                mobile_geo_name = str(new_mobile_obj.subunits) + ',' + str([int(x) for x in mobile_geo])
                number_neighbors[stat_geo_name].add(mobile_geo_name)
                if mobile_geo_name not in number_neighbors.keys():
                    number_neighbors[mobile_geo_name] = set()
                    number_neighbors[mobile_geo_name].add(stat_geo_name)
                
                #assoc is used to remind the clash check that static objects should be removed from the grid before
                #deleting any of the mobile objects.
                associated_statics.append((mobile_geo_name, copy_transform_object(static_object)))

                append_max -= 1
                append_obj += 1
                if static_main == 'A_D':
                    break
                else:
                    print static_axis_used
                    if layer != 1:
                        if static_axes[0] == compatible_axes[static_axis_used][-1]:
                            break
            '''
            if layer < layers:
                #always add the number of appended objects to the branch_sum
               branch_sum += float(append_obj)
            '''
        #perform clash check with the final_scaffold_obj and objects in last_layer_placed
        print new_mobile_obj
        if not last_layer_placed:
            break
        pop_count = 0
        #this pre_count is set up to figure out just how large the accepted list is so far
        #it should be able to be done away with and just use the len(obj_in_layer)
        #pre_count_of_accepted = len(list_of_accepted_additions)
        grid_for_current_layer = set()
        for index, obj_in_layer in enumerate(last_layer_placed[:]):
            #obj_in_layer.write_to_file(str(layer) + '_' + str(index)+'.pdb')
            #in associ. 0 is the name in connection dict and 1 is the object
            copy_static = copy_transform_object(associated_statics[index][1])# - pop_count][1])
            copy_static_grid = copy_static.get_atom_grid(grid_size=clash_mag)
            #take out the static object that the mobile is placed next to when checking if it should stay
            working_grid = final_scaffold_grid - copy_static_grid

            copy_obj = copy_transform_object(obj_in_layer[0])
            copy_obj_grid = copy_obj.get_atom_grid(grid_size=clash_mag)
            grid_check = copy_obj_grid & working_grid
            if grid_check: #if grid_check has clashing members
                #last_layer_placed.pop(index - pop_count)
                try:
                    del number_neighbors[associated_statics[index]]# - pop_count][0]]
                except KeyError:
                    pass
                for k, v in number_neighbors.iteritems():
                    v.discard(associated_statics[index])# - pop_count][0])

                kill_index = last_layer_placed.index(obj_in_layer)
                del last_layer_placed[kill_index]
                #del last_layer_placed[index - pop_count]
                #pop_count += 1
                print 'clash'
            else:
                final_scaffold_grid = final_scaffold_grid | copy_obj_grid
                #list_of_accepted_additions.append(copy_transform_object(copy_obj))
                copy_obj.append_to_file(temp_file_name, model=True)
                del copy_obj

                distances.append(rigid_body_transformation[index][0])
                rotations.append(rigid_body_transformation[index][1])
                print 'non_clash'

                grid_for_current_layer = grid_for_current_layer | copy_obj_grid
        #layer 0 and all evens are A additions
        #layer 1 and all odss are C additions
        if layer % 2 == 0:
            a_count += len(last_layer_placed) #list_of_accepted_additions) - pre_count_of_accepted
            if layer < layers:
                a_not_last_layer += len(last_layer_placed) #list_of_accepted_additions) - pre_count_of_accepted
        else:
            b_count += len(last_layer_placed) #list_of_accepted_additions) - pre_count_of_accepted
            branch_sum += len(last_layer_placed) #list_of_accepted_additions) - pre_count_of_accepted


        if particle_goal < a_count + b_count:
            break
    #convert final grid to xyz coords to get relative size of assembly
    #print len(grid_for_current_layer)
    #grid_diameter = find_longest_grid_diameter(grid_for_current_layer)
    grid_diameter = 0
    ratio = float(branch_sum)/float(a_not_last_layer)
    csv_line = file_name + ',' + str(layer) + ',' + str(init_null) + ',' + str(init_frac) +',' +\
               str(a_count) + ',' + str(b_count) + ',' + str(ratio) + ',' + str(int(grid_diameter)) + '\n'
    new_file_name = file_name + '_layer' + str(layer) + '_null' + str(init_null) +\
                            '_frac' + str(init_frac) + '_A' + str(a_count) + '_B' + str(b_count) +\
                            '_R' + str(ratio) + '_size' + str(int(grid_diameter)) + '.pdb'
    #for obj in list_of_accepted_additions:
    #obj.append_to_file(file_name, model=True)
    os.rename(temp_file_name, new_file_name)

    #del list_of_accepted_additions

    with open('all_fractal_data_null%s_frac%s.csv' % (init_null, init_frac), 'a') as mycsv:
        mycsv.write(csv_line)
    with open('distances_%s.txt' % new_file_name[:-4], 'w') as mydist:
        for line in distances:
            mydist.write(str(line) + '\n')
    with open('rotations_%s.txt' % new_file_name[:-4], 'w') as myrots:
        for line in rotations:
            myrots.write(str(line) + '\n')
    with open('neighbors_%s.txt' % new_file_name[:-4], 'w') as myneigh:
        for k, v in number_neighbors.iteritems():
            myneigh.write(k + ',' + str(len(v)) + '\n')

   # for name, my_obj in collections.OrderedDict(sorted(Two_comp_obj.ensemble.items())).iteritems():
   #     if name.split('_')[-1].startswith('trans'):
   #         my_obj.write_to_file('%s.pdb' % name, 'A')
    del final_scaffold_grid 
    
    return

################################################################################
#
#
################################################################################

def main(argv):
    
    parser = optparse.OptionParser(usage="\n\nGenerate csv of symmetric cofactor matches.")
    
    parser.add_option('--static', dest = 'static_filename',
        help = 'The static object')
    
    parser.add_option('--mobile', dest = 'mobile_filename',
        help = 'The mobile objcet which needs to be d2 for now')
    
    parser.add_option('--layers', type="int", nargs=1, dest = 'layers',
        default = 1,
        help = 'Number of layers appended to static seed ( default = 1 )')
    
    parser.add_option('--mobile-axis-bin', type="str", dest = 'mobile_axis_bin',
        default = [],
        help = 'Comma separated bin of possible c-symm axes for mobile object( default = '' )')
    
    parser.add_option('--static-axis-bin', type="str", dest = 'static_axis_bin',
        default = [],
        help = 'Comma separated bin of possible c-symm axes for static object( default = '' )')
    
    parser.add_option('--trans-bin', type="str", dest = 'trans_bin',
        default = [],
        help = 'Comma separated bin of possible translations along c-symm axis( default = '' )')
    
    parser.add_option('--trans-prob', type="str", dest = 'trans_prob',
        default = [],
        help = 'Comma separated set of probabilities corresponding to trans-bin ( default = '' )')
    
    parser.add_option('--rot-bin', type="str", dest = 'rot_bin',
        default = [],
        help = 'Comma sepparated bin of possible rotations about c-symm axis ( default = '' )')

    parser.add_option('--rot-prob', type="str", dest = 'rot_prob',
        default = [],
        help = 'Comma separated set of probabilities corresponding to rot-bin ( default = '' )')

    parser.add_option('--null1', type="float", nargs=1, dest = 'null1',
        default = 0.0,
        help = 'Probability of not appending first scaffold, must be <1.0 ( default = 0.0 )')

    parser.add_option('--null2', type="float", nargs=1, dest = 'null2',
        default = 0.0,
        help = 'Probability of not appending second scaffold, must be <1.0 ( default = 0.0 )')
    
    parser.add_option('--paired-bins', action="store_true", dest = 'paired_bins',
        default = False,
        help = 'Translation and rotation values are paired ( default = 1 )')
    
    parser.add_option('--paired-prob', type="str", dest = 'paired_prob',
        default = [],
        help = 'Comma separated set of probabilities corresponding to paired trans and rot bins ( default = '' )')
    
    parser.add_option('--clash-mag', type="float", nargs=1, dest = 'clash_mag',
        default = 3.0,
        help = 'Minimum atom-atom clash distance ( default = 3.0 )')
    
    parser.add_option('--a-frac', type="float", nargs=1, dest = 'a_frac',
        default = 0.0,
        help = 'Fraction of first scaffold selected for each layer, must be <1.0 ( default = 0.0 )')
    
    parser.add_option('--b-frac', type="float", nargs=1, dest = 'b_frac',
        default = 0.0,
        help = 'Fraction of first scaffold selected for each layer, must be <1.0 ( default = 0.0 )')
    
    parser.add_option('--out-file-name', dest = 'file_name',
        help = 'The desired file name output.')
    
    parser.add_option('--D2-main-axis', type="str", dest = 'd2_main_axis',
        default = [],
        help = 'Sets a particular axis for any D2 proteins as main axis( default = '' )')
    
    (options,args) = parser.parse_args()
    
    #time1 = time.clock() 
    static = read_file(str(options.mobile_filename))
    mobile = read_file(str(options.static_filename))
    

    if options.trans_bin:
        trans_bin = options.trans_bin.split(',')# [float(x) for x in options.trans_bin.split(',')]
        trans_range = trans_bin[0].split('-')
        if len(trans_range) == 2:
            trans_bin = [ float(x) for x in range(int(float(trans_range[0])),\
                                                  int(float(trans_range[1])),\
                                                  int(float(trans_bin[1]))) ] 
        else:
            trans_bin = [float(x) for x in trans_bin]
    else:
        trans_bin = []
    if options.trans_prob:
        trans_prob = [float(x) for x in options.trans_prob.split(',')]
    else:
        trans_prob = []
    if options.rot_bin:
        rot_bin = options.rot_bin.split(',') # [float(x) for x in options.rot_bin.split(',')]
        rot_range = rot_bin[0].split('-')
        if len(rot_range) == 2:
            rot_bin = [ float(x) for x in  range(int(float(rot_range[0])), \
                                                 int(float(rot_range[1])), \
                                                 int(float(rot_bin[1]))) ]
        else:
            rot_bin = [float(x) for x in rot_bin]
    else:
        rot_bin = []
    if options.rot_prob:
        rot_prob = [float(x) for x in options.rot_prob.split(',')]
    else:
        rot_prob = []
    
    if not options.paired_bins:
        options.paired_bins = False

    if options.paired_prob:
        paired_prob = [float(x) for x in options.paired_prob.split(',')]
    
    if options.static_axis_bin:
        static_axis_bin = options.static_axis_bin.split(',')
    if options.mobile_axis_bin:
        mobile_axis_bin = options.mobile_axis_bin.split(',')
    '''
    dn_dn_layerer(options.static_filename, options.mobile_filename, options.layers,\
                  trans_bin, trans_prob, rot_bin, rot_prob, options.null1, options.null2,\
                  options.paired_bins, paired_prob, options.d2_main_axis,\
                  options.clash_mag, options.a_frac, options.b_frac,\
                  options.file_name, static_axis_bin, mobile_axis_bin)
    '''
    #null_set = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    #frac_set = [1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
    null_set = [0.0]#,0.1,0.2,0.3,0.4]
    frac_set = [1.0]#,0.9,0.8,0.7,0.6]
    for null in null_set:
        for frac in frac_set:
            dn_dn_layerer(options.static_filename, options.mobile_filename, options.layers,\
                          trans_bin, trans_prob, rot_bin, rot_prob, options.null1, null,\
                          options.paired_bins, paired_prob, options.d2_main_axis,\
                          options.clash_mag, options.a_frac, frac,\
                          options.file_name, static_axis_bin, mobile_axis_bin)
            print "We are on:", null, frac 
    
    #print 'It took %.4f seconds to run the Dn_Dn_assem_maker program.' % float(time.clock() - time1)
    
if __name__ == "__main__":
    main(sys.argv[1:])
