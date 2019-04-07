#!/usr/bin/env python
from transform import *
import math, sys, os, time, random

class TwoComponentAssem(object):
    
    def __init__(self, static_pdb_path, mobile_pdb_path):
        if type(static_pdb_path) == str:
            self.static_pdb = read_file(static_pdb_path)
            self.static_obj = Transform(self.static_pdb)
        else:
            print "TwoComponentAssem arguement 1 is not of type str\n \
                   and thus not a valid path name."
            sys.exit()

        if type(mobile_pdb_path) == str:
            self.mobile_pdb = read_file(mobile_pdb_path)
            self.mobile_obj = Transform(self.mobile_pdb)
        else:
            print "TwoComponentAssem arguement 2 is not of type str\n \
                   and thus not a valid path name"
            sys.exit()
        self.alignment_dictionary = {}
        self.ensemble = {}
        self.static_clash_chains = set()
        self.mobile_clash_chains = set()
    
    def swap_object_properties(self):
        new_static = copy_transform_object(self.mobile_obj)
        new_mobile = copy_transform_object(self.static_obj)
        self.static_obj = new_static
        self.mobile_obj = new_mobile

    def empty_alignment_dictionary(self):
        self.alignment_dictionary = {}

    def empty_ensemble(self):
        self.ensemble = {}

    def get_axis_coord_by_name(self, name, static=False):
        if static == False:
            axes = self.mobile_obj.get_symmetric_axes()
            return axes[name]
        else:
            axes = self.static_obj.get_symmetric_axes()
            return axes[name]
    
    def get_vector_difference(self):
        return self.mobile_obj.get_geo_center() - self.static_obj.get_geo_center()

    def origin_align(self, flip=False, axis1='', axis2=''):
        self.static_obj.translate([0.,0.,0.])
        self.mobile_obj.translate([0.,0.,0.])
        #static_axes = self.static_obj.get_symmetric_axes()
        #mobile_axes = self.mobile_obj.get_symmetric_axes()
        self.mobile_obj.major_symmetric_axis_align(self.static_obj, flip)#, axis1, axis2)
        #else:
        #    self.mobile_obj.major_symmetric_axis_align(self.static_obj, flip)

    def align_chain_A(self, inverse=False):
        static_axes = self.static_obj.get_symmetric_axes()
        mobile_axes = self.mobile_obj.get_symmetric_axes()
        self.mobile_obj.align_object_to_vector(mobile_axes['A'], static_axes['A'], True)
    
    # Can pass an esemble to this definition, if ensemble is a 
    # dictionary it will unpack the dictionary. If true it defaults
    # to whichever is true first (self.ensemble > self.alignment_dictionary)
    # if left as false it will translate the mobile object along axis
    def translate_along_axis(self, axis, interval=1.0, min_distance=0.0,\
                             max_distance=0, ensemble=False,\
                             check_clash=False, output_name=''):
        ###FINISH THIS
        if int(interval) == int(min_distance) == int(max_distance) == 0:
            return
        mobile_objects = {}
        axis_to_origin = get_vector_difference(axis, self.mobile_obj.get_geo_center()[:])
        axis_mag = get_mag(axis_to_origin)
        unit_vector = [x/axis_mag for x in axis_to_origin]
        if ensemble == False:
            if output_name:
                mobile_objects[output_name] = copy_transform_object(self.mobile_obj)
            else:    
                mobile_objects['mobile_obj'] = copy_transform_object(self.mobile_obj)
        else:
            if type(ensemble) == dict:
                operating_ensemble = ensemble
            elif self.ensemble:
                operating_ensemble = self.ensemble
            elif self.alignment_dictionary:
                operating_ensemble = self.alignment_dictionary
            else:
                print "You marked ensemble True in 'translate_along_axis' yet\
                       self.ensemble and self.alignment_dictionary are empty"
            for name, ensemble_obj in operating_ensemble.iteritems():
                mobile_objects[name] = copy_transform_object(ensemble_obj)
        static_copy = copy_transform_object(self.static_obj)
        static_copy.get_xyz_by_atomtype(['N','CA','C','O'], True)
        '''if self.static_clash_chains:
            static_copy.get_xyz_by_chain(list(self.static_clash_chains), True)
            static_copy.get_xyz_by_atomtype(['N','CA','C','O'], True)
        else:
            static_copy.get_xyz_by_atomtype(['N','CA','C','O'], True)'''
        
        mobile_center = self.mobile_obj.get_geo_center()
        for name, copy_obj in mobile_objects.iteritems():
            translation = min_distance
            # ensembles will have more than 1 dictionary item so make copyies
            if max_distance == 0:
                while copy_obj.clash_check(static_copy, 3.0, ['N','CA','C','O'],\
                                      [], list(self.mobile_clash_chains))\
                                      == True:
                    copy_obj.translate(add_vector_difference([x*translation for x in unit_vector],\
                                       self.mobile_obj.get_geo_center()[:]))
                    translation += interval
                self.ensemble['%s_trans%.4f' % (name, translation)] \
                                            = copy_transform_object(copy_obj)
            else:
                print 'translation'
                while translation <= max_distance:
                    print 'trans', translation
                    print 'u-vec', unit_vector
                    copy_obj.translate(add_vector_difference([x*translation for x in unit_vector],\
                                       mobile_center[:]))
                    if check_clash:
                        if copy_obj.clash_check(static_copy, 3.0, ['N','CA','C','O'],\
                                           [], []) == False:
                            self.ensemble['%s_trans%.4f' % (name, translation)] \
                                                = copy_transform_object(copy_obj)
                    else:
                        self.ensemble['%s_trans%.4f' % (name, translation)] \
                                            = copy_transform_object(copy_obj)
                        if len(mobile_objects.keys()) == 1:
                            self.mobile_obj = copy_transform_object(copy_obj)
                    translation += interval
                    
        
    # Can pass an esemble to this definition, if ensemble is a 
    # dictionary it will unpack the dictionary. If true it defaults
    # to whichever is true first (self.ensemble > self.alignment_dictionary)
    # it is advised to either provide the alignment ensemble or empty
    # ensemble before using this definition for rotating about obj's own axis
    # if left as false it will translate the mobile object along axis
    def rotate_about_axis(self, axis, interval=10.0, min_rotation=0.0,\
                             max_rotation=360.0, ensemble=False,\
                             check_clash=False, output_name=''):
        ###FINISH THIS
        mobile_objects = {}
        if ensemble == False:
            if output_name:
                mobile_objects[output_name] = copy_transform_object(self.mobile_obj)
            else:    
                mobile_objects['mobile_obj'] = copy_transform_object(self.mobile_obj)
        else:
            if type(ensemble) == object:
                operating_ensemble = ensemble
            elif self.ensemble:
                operating_ensemble = self.ensemble
            elif self.alignment_dictionary:
                operating_ensemble = self.alignment_dictionary
            else:
                print "You marked ensemble True in 'rotate_about_axis' yet\
                       self.ensemble and self.alignment_dictionary are empty"
            for name, ensemble_obj in operating_ensemble.iteritems():
                mobile_objects[name] = copy_transform_object(ensemble_obj)

        if check_clash == True:
            static_copy = copy_transform_object(self.static_obj)
            if self.static_clash_chains:
                static_copy.get_xyz_by_chain(list(self.static_clash_chains), True)
                static_copy.get_xyz_by_atomtype(['N','CA','C','O'], True)
        
        for name, copy_obj in mobile_objects.iteritems():
            rotation = min_rotation
            # ensembles will have more than 1 dictionary item so make copyies
            print 'rotation'
            ref_copy_obj = copy_transform_object(copy_obj)
            while rotation < max_rotation:
                copy_obj.rotate_object(list(copy_obj.get_geo_center()),\
                                       list(axis), math.radians(rotation))
                if check_clash == True:
                    if copy_obj.clash_check(static_copy, 3.0, ['N','CA','C','O'],\
                                            [], list(self.mobile_clash_chains))\
                                            == False:
                        self.ensemble['%s_rot%.4f' % (name, rotation)] \
                                            = copy_transform_object(copy_obj)
                else:
                    self.ensemble['%s_rot%.4f' % (name, rotation)] \
                                            = copy_transform_object(copy_obj)
                        
                copy_obj = copy_transform_object(ref_copy_obj)
                rotation += interval

    #This assumes that the static and mobile object are in the correct binding orientation for one chain
    def add_asymmetric_mobile_to_static(self, null=0.0):
        original_static_location = self.static_obj.get_geo_center()[:]
        reverse_static_loc = [-1.*x for x in original_static_location]
        #origin align here
        ref_static = copy_transform_object(self.static_obj)
        self.static_obj.translate([0.,0.,0.], update=True)
        
        ref_mobile = copy_transform_object(self.mobile_obj)
        
        #translate the mobile the same distance and direction of the static
        for index, xyz in enumerate(self.mobile_obj.xyz_total):
            self.mobile_obj.xyz_total[index] = add_vector_difference(xyz[:], reverse_static_loc[:])

        symm_static_matrix = self.static_obj.calculate_symmetric_transform()
        transformed_objects = self.mobile_obj.apply_transform_to_object(symm_static_matrix)
        random.shuffle(transformed_objects)
        

        final_objs = []
        #pop_count = 0
        for ind, obj in enumerate(transformed_objects):
            #pick a value from trans_bin and rot_bin, then apply to mobile.
            if random.random() < null:
                #transformed_objects.pop(ind - pop_count)
                #pop_count += 1
                continue
            else:
                for index, xyz in enumerate(obj.xyz_total):
                    obj.xyz_total[index] = add_vector_difference(xyz[:], original_static_location)
                final_objs.append(copy_transform_object(obj))
        
        self.mobile_obj = copy_transform_object(ref_mobile)
        self.static_obj = copy_transform_object(ref_static)
        
        return final_objs

class Dn_D2_AssemMaker(TwoComponentAssem):
    
    def __init__(self, static_pdb_path, mobile_pdb_path):
        if type(static_pdb_path) == str and type(mobile_pdb_path) == str:
            self.static_pdb_path = static_pdb_path
            self.mobile_pdb_path = mobile_pdb_path
            self.alignment_dictionary = {}
            self.ensemble = {}
            self.static_clash_chains = set()
            self.mobile_clash_chains = set()
            super(Dn_D2_AssemMaker, self).__init__(self.static_pdb_path, self.mobile_pdb_path)
        
        else:
            print "Dn_D2_AssemMaker arguement 1 or 2 is not of type str\n \
                   and thus not a valid path name."
            sys.exit()
    
    def get_dn_c2_axes_names(self, static=False):
        if static == True:
            axes_dict = self.static_obj.get_symmetric_axes()
        else:
            axes_dict = self.mobile_obj.get_symmetric_axes()
        trans_axes_names = []
        
        for key in axes_dict.keys():
            print key
            if len(key.split('_')) == 2: #and key[1] == '_':
                '''chains = key.split('_')
                if 'A' in chains or '1' in chains:
                    if static == True:
                        self.static_clash_chains.add(chains[0])
                        self.static_clash_chains.add(chains[1])
                    else:
                        self.mobile_clash_chains.add(chains[0])
                        self.mobile_clash_chains.add(chains[1])'''
                val = axes_dict[key]
                trans_axes_names.append(str(key))
        print 'important axes names: ', trans_axes_names
        return trans_axes_names

    def generate_alignment_dictionary(self, static_axes=[], mobile_axes=[], invert=False, static_main_axis='', mobile_main_axis=''):
        #The idea here is to make copies and fill alignment_dictionary
        #First determine which is the D2 protein
        print 'Starting generate_alignment_dictionary'
        #self.static_obj.translate([0.,0.,0.])
        #self.mobile_obj.translate([0.,0.,0.])
        
        self.origin_align(flip=False, axis1=mobile_main_axis, axis2=static_main_axis)
        static_axes_coordinates = self.static_obj.get_symmetric_axes()
        mobile_axes_coordinates = self.mobile_obj.get_symmetric_axes()
        static_c2_AxesOI = self.get_dn_c2_axes_names(static=True)
        mobile_c2_AxesOI = self.get_dn_c2_axes_names(static=False)
        if static_axes:
            static_c2_AxesOI = static_axes
        if mobile_axes:
            mobile_c2_AxesOI = mobile_axes
        for static_axis in static_c2_AxesOI:
            for mobile_axis in mobile_c2_AxesOI:
                ref_copy = copy_transform_object(self.mobile_obj)
                try:
                    if invert == False:
                        ref_copy.align_object_to_vector(axis1=mobile_axes_coordinates[mobile_axis],\
                                                        axis2=static_axes_coordinates[static_axis],\
                                                        update=True)
                    else:
                        ref_copy.align_object_to_vector(axis1=[-1.*x for x in mobile_axes_coordinates[mobile_axis]],\
                                                        axis2=static_axes_coordinates[static_axis],\
                                                        update=True)
                except KeyError:
                    print 'Key error in generate alignment ensemble'
                    return False
                '''mobile_axes_coordinates = self.mobile_obj.get_symmetric_axes()
                ortho_planar_vector = cross(mobile_axes_coordinates['minor'], \
                                            mobile_axes_coordinates[mobile_axis])
                self.mobile_obj.align_object_to_vector(axis1=ortho_planar_vector,\
                                                       axis2=static_axes_coordinates['major'],\
                                                       update=True)'''
                self.alignment_dictionary['%s_aligned_to_%s' %\
                                      (mobile_axis, static_axis)] = copy_transform_object(ref_copy)
                #self.mobile_obj = copy_transform_object(ref_copy)
                return True
    
    def generate_c2_axes_positional_ensemble(self, rotation_min=0, rotation_max=360.0,\
                                             rotation_interval=10.0, translation_min=0.0,\
                                             translation_max=0.0, translation_interval=1.0,\
                                             static_axes=[], mobile_axes=[], \
                                             static_alignment_axis='', mobile_alignment_axis=''):
        master_ensemble = {}
        self.generate_alignment_dictionary(static_axes, mobile_axes, False, static_alignment_axis, mobile_alignment_axis)
        ref_obj = copy_transform_object(self.mobile_obj)
        for name, obj in self.alignment_dictionary.iteritems():
#Below the code is fine, just set it up properly for dn_dn with the if statement
            if self.mobile_obj.subunits in [4,6]:
                axes_aligned = name.split('_aligned_to_')
                d2_obj_axes_coords = obj.get_symmetric_axes()
                self.mobile_obj = copy_transform_object(obj)
                self.rotate_about_axis(d2_obj_axes_coords[axes_aligned[0]],\
                                       rotation_interval, rotation_min, rotation_max,\
                                       False, False, name) #ensemble and check_clash are false
                #check here
                self.translate_along_axis(d2_obj_axes_coords[axes_aligned[0]],\
                                          translation_interval, translation_min,\
                                          translation_max, self.ensemble, False) #made this clash free
                print 'came out of while loop'
                for new_name, final_obj in self.ensemble.iteritems():
                    master_ensemble[new_name] = final_obj
                self.ensemble = {}
            self.mobile_obj = copy_transform_object(ref_obj)
        self.ensemble = master_ensemble
        #print master_ensemble.keys()
        #print self.ensemble.keys()
        #sys.exit()

    #null must be a number between 0.0 and 1.0 and is the probability that nothing will be placed for a given iteration
    #translation_bin are the possible translations along a csymmetric axis one can make
    #rotation_bin are the possible rotations about the csymmetric axis one can make
    #probability lists correspond to the weight of each translation/rotation otherwise each value you will be equal
    def add_mobile_to_static(self, translation, rotation, null=0.0,\
                             static_axis='', mobile_axis='',\
                             static_main='', mobile_main=''):
        
        if null> 0.0:
            if random.random() < null:
                print 'nulled out'
                return []
        
        #make a copy of mobile and locate all termini
        ref_obj = copy_transform_object(self.mobile_obj)
        original_static_location = self.static_obj.get_geo_center()[:]
        #origin align here
        self.alignment_dictionary = {}
        check_key = self.generate_alignment_dictionary([static_axis], [mobile_axis], invert=True, \
                                           mobile_main_axis=mobile_main, static_main_axis=static_main)
        if check_key == False:
            print 'Found a key error, in generate_alignment'
            return None

        self.mobile_obj = copy_transform_object(self.alignment_dictionary[mobile_axis + '_aligned_to_' + static_axis])
        d2_obj_axes_coords = self.mobile_obj.get_symmetric_axes()
        print d2_obj_axes_coords
        #self.mobile_obj = copy_transform_object(obj)
        mobile_center = self.mobile_obj.get_geo_center()
        try:
            self.mobile_obj.rotate_object(mobile_center[:], \
                                          d2_obj_axes_coords[mobile_axis], \
                                          math.radians(rotation))
            self.mobile_obj.translate_along_axis(d2_obj_axes_coords[mobile_axis], \
                                                 magnitude=translation, invert=True)
        except KeyError:
            print 'Found a key error, in d2_obj_axes_coords'
            self.mobile_obj = copy_transform_object(ref_obj)
            return None
        for index, xyz in enumerate(self.mobile_obj.xyz_total):
            self.mobile_obj.xyz_total[index] = add_vector_difference(original_static_location, xyz) 

        output_mobile = copy_transform_object(self.mobile_obj)
        self.mobile_obj = copy_transform_object(ref_obj)
        

        return output_mobile
    
    def add_mobile_to_static_bu(self, trans_bin=[], trans_prob=[],  rot_bin=[], rot_prob=[],\
                             null=0.0, paired_bins=False, static_axes=[], mobile_axes=[],\
                             static_main='', mobile_main='', paired_prob=''):
        #make a copy of mobile and locate all termini
        original_static_location = self.static_obj.get_geo_center()[:]
        #origin align here
        self.alignment_dictionary = {}
        static_set = set()
        for axes in static_axes:
            for axis in axes:
                static_set.add(axis)
        static_list = list(static_set)
        mobile_set = set()
        for axes in mobile_axes:
            for axis in axes:
                mobile_set.add(axis)
        mobile_list = list(mobile_set)

        #right here it is passing an empty list for one or more of the axes
        self.generate_alignment_dictionary(static_list, mobile_list, invert=True, \
                                           mobile_main_axis=mobile_main, static_main_axis=static_main)
        ref_obj = copy_transform_object(self.mobile_obj)
        random.shuffle(static_list)
        output_objects = []
        if paired_bins == True:
            axis_trans_rot_prob = zip(static_axes, mobile_axes, trans_bin, rot_bin, paired_prob)
            rand_frac = random.random()
            tuple_probs = zip(paired_prob,range(len(paired_prob)))
            tuple_probs.sort()
            sum_of_probs = 0.0
            enumerated_probs = []
            for prob, ind in tuple_probs:
                if prob + sum_of_probs > rand_frac:
                    translation = trans_bin[ind]
                    rand_index = ind
                    break
                    sum_of_probs += prob
            
        
        for axis_type in static_list:
            
            if axis_type == terminating_axis:
                if non_term_placed:
                    continue
                    
            #for name, obj in self.alignment_dictionary.iteritems():
            #pick a value from trans_bin and rot_bin, then apply to mobile.
            #if random.random() < null:
                #self.alignment_dictionary[name] = []
                #self.mobile_obj = copy_transform_object(ref_obj)
                #continue
            if not trans_prob:
                rand_index = random.randrange(0,len(trans_bin))
                print rand_index
                print len(trans_bin)
                translation = trans_bin[rand_index]
            else:
                rand_frac = random.random()
                tuple_probs = zip(trans_prob,range(len(trans_prob)))
                tuple_probs.sort()
                sum_of_probs = 0.0
                enumerated_probs = []
                for prob, ind in tuple_probs:
                    if prob + sum_of_probs > rand_frac:
                        translation = trans_bin[ind]
                        rand_index = ind
                        break
                    sum_of_probs += prob

            if paired_bins == True:
                rotation = rot_bin[rand_index]
                check_axes = name.split('_aligned_to_')
                if check_axes[0] not in mobile_axes[rand_index]:
                    continue
                if check_axes[1] not in static_axes[rand_index]:
                    continue
            else:
                if random.random() < null:
                    continue
                if not rot_prob:
                    rand_rot_index = random.randrange(0,len(rot_bin))
                    rotation = rot_bin[rand_rot_index]
                else:
                    rand_frac = random.random()
                    tuple_probs = zip(rot_prob,range(len(rot_prob)))
                    tuple_probs.sort()
                    sum_of_probs = 0.0
                    enumerated_probs = []
                    for prob, ind in rot_probs:
                        if prob + sum_of_probs > rand_frac:
                            rotation = rot_bin[ind]
                            break
                        sum_of_probs += prob
            #Below the code is fine, just set it up properly for dn_dn with the if statement
            #if self.mobile_obj.subunits == 4:
            axes_aligned = name.split('_aligned_to_')
            print 'axes_aligned:', axes_aligned
            d2_obj_axes_coords = obj.get_symmetric_axes()
            print d2_obj_axes_coords
            self.mobile_obj = copy_transform_object(obj)
            try:
                self.mobile_obj.rotate_object(self.mobile_obj.get_geo_center()[:], \
                                              d2_obj_axes_coords[axes_aligned[0]], \
                                              math.radians(rotation))
                self.mobile_obj.translate_along_axis(d2_obj_axes_coords[axes_aligned[0]], \
                                                     magnitude=translation, invert=True)
            except KeyError:
                self.mobile_obj = copy_transform_object(ref_obj)
                continue
            for index, xyz in enumerate(self.mobile_obj.xyz_total):
                self.mobile_obj.xyz_total[index] = add_vector_difference(original_static_location, xyz) 

            output_objects.append(copy_transform_object(self.mobile_obj))
            self.mobile_obj = copy_transform_object(ref_obj)
        
            if name == terminating_axis:
                break
            else:
                non_term_placed = True

        return output_objects
                                
class Cn_Cn_AssemMaker(TwoComponentAssem):
    
    def __init__(self, static_pdb_path, mobile_pdb_path):
        if type(static_pdb_path) == str and type(mobile_pdb_path) == str:
            self.static_pdb_path = static_pdb_path
            self.mobile_pdb_path = mobile_pdb_path
            self.alignment_dictionary = {}
            self.ensemble = {}
            self.static_clash_chains = set()
            self.mobile_clash_chains = set()
            super(Cn_Cn_AssemMaker, self).__init__(self.static_pdb_path, self.mobile_pdb_path)
        
        else:
            print "Cn_Cn_AssemMaker arguement 1 or 2 is not of type str\n \
                   and thus not a valid path name."
            sys.exit()
    
    def generate_alignments(self):
        #The idea here is to make copies and fill alignment_dictionary
        #First determine which is the D2 protein
        self.origin_align()
        self.align_chain_A()
        self.alignment_dictionary['M1M1_A'] = copy_transform_object(self.mobile_obj)
        self.origin_align(True)
        self.align_chain_A()
        self.alignment_dictionary['M1M2_A'] = copy_transform_object(self.mobile_obj)

    def generate_conformation(self, T=0.0, IR=0.0, ER=0.0, O=0.0, check_clash=False, origin_align=False, invert_mobile=False, invert_translation=False, return_key=False):
        reference_obj = copy_transform_object(self.mobile_obj)
        
        if origin_align == True:
            self.origin_align()
        if invert_mobile == True:
            #self.mobile_obj.rotate_object([0.0,0.0,0.0], self.mobile_obj.get_symmetric_axes()['major'], math.pi)
            self.mobile_obj.rotate_object([0.0,0.0,0.0], self.mobile_obj.major_symm_axis, math.pi)
            
        self.mobile_obj.translate_along_axis('A', T, invert=invert_translation)

        if IR != 0.0:
            mobile_obj_axes_coords = self.mobile_obj.get_symmetric_axes()
            #self.mobile_obj.rotate_object(self.mobile_obj.get_geo_center(), mobile_obj_axes_coords['major'], IR, True)
            self.mobile_obj.rotate_object(self.mobile_obj.get_geo_center(), self.mobile_obj.major_symm_axis, IR, True)
            
        if ER != 0.0:
            static_obj_axes_coords = self.static_obj.get_symmetric_axes()
            #self.mobile_obj.rotate_about_axis(static_obj_axes_coords['major'], ER, True)
            self.mobile_obj.rotate_about_axis(self.static_obj.major_symm_axis, ER, True)
        
        #self.mobile_obj.translate_along_axis('major', O)
        self.mobile_obj.translate_along_axis(self.static_obj.major_symm_axis, O)

        ensemble_name = 'T%.2f_IR%.2f_ER%.2f_O%.2f' % (T,IR,ER,O)

        if check_clash == True:
            if not self.mobile_obj.clash_check(self.static_obj):
                self.ensemble[ensemble_name] = copy_transform_object(self.mobile_obj)
                self.mobile_obj = copy_transform_object(reference_obj)
                if return_key == True:
                    return ensemble_name
        else:
            self.ensemble[ensemble_name] = copy_transform_object(self.mobile_obj)
            self.mobile_obj = copy_transform_object(reference_obj)
            if return_key == True:
                return ensemble_name
