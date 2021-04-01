#!/usr/bin/env python
# coding: utf-8

# In[1]:



vcf_file_address = "22_5sample.vcf"  # "22_4k.vcf" 
vcf_file = open(vcf_file_address,'r');


blockid_blocks_s1 = []
allele_blocks_s1 = []
varpos_blocks_s1 = []


blockid_blocks_s2 = []
allele_blocks_s2 = []
varpos_blocks_s2 = []

blockid_blocks_s3 = []
allele_blocks_s3 = []
varpos_blocks_s3 = []

blockid_blocks_s4 = []
allele_blocks_s4 = []
varpos_blocks_s4 = []

var_blockids1_blocks_s2 = []   # for each var in s2, shows the s1 blockid
var_blockids1_blocks_s3 = []   # for each var in s3, shows the s1 blockid
var_blockids1_blocks_s4 = []   # for each var in s4, shows the s1 blockid


# for parental that it's only one block per chr so it doesn't have blockid
alleles_s5 = []
varposs_s5 = []


varposs_all = []

num_weird_occurs=0

var_lines_list = []

header_lines_list = []
for line in vcf_file:
    line_strip = line.strip()
    
    if line_strip.startswith('#'):
        header_lines_list.append(line_strip)
        sample_names = line_strip.split('\t')[9:]   # last line of header contains sample name
        #print(sample_names)
        
    else:
        line_parts = line_strip.split('\t')
        var_lines_list.append(line_parts)
        #gt_flags, s1, sss, ss, s2, s5 = line_parts[8:]
        gt_flags, s1, s2,s3,s4, s5 = line_parts[8:]       # hifi-whatshp illumina-whatshp hifi-hapcut illumina-hapcut true

        varpos = int(line_parts[1])                 # variant position
        
        varposs_all.append(varpos)
        gt_flags_split = gt_flags.split(":")
        
        if "PS" in gt_flags_split:
            idx_ps = gt_flags.split(":").index("PS")
            idx_gt = gt_flags.split(":").index("GT")
            

                                                 ## s5- parental - last sample in the vcf file
            split_s5 = s5.split(":")
            allele_s5 = split_s5[idx_gt]
            blockid_s5 = split_s5[idx_ps]
            
            if allele_s5 == "0|1" or allele_s5 == "1|0":
                alleles_s5.append(int(allele_s5[0]))
                varposs_s5.append(varpos)

            
            
                                                ###  sample 1 ###
            split_s1 = s1.split(":")
            blockid_s1 = split_s1[idx_ps]
            allele_s1 = split_s1[idx_gt]
            

            if blockid_s1 == ".":                                 # unphased variants for sample one
                pass
            
            elif blockid_s1 in blockid_blocks_s1:
                assert allele_s1[1] == "|", "allele_s1[1] should be | and so phased. A non-phased variant is in the phase block"
                
                #if blockid_s1 == previous_blockid_s1: previous_blockid_s1_idx has been defined in previous variant which is the last idx
                if blockid_s1 != previous_blockid_s1:  
                    # weird situation happens for whatshap in which two variants are in the same block but in between there is/are other phaseblock(s)
                    # note that the unphased variants are ignored sooner and are not considered here
                    previous_blockid_s1_idx = blockid_blocks_s1.index(blockid_s1)
                    #print("two variants are in the same block, but far from each other",varpos)                
                    num_weird_occurs+=1

                allele_blocks_s1[previous_blockid_s1_idx].append( int(allele_s1[0]) )
                varpos_blocks_s1[previous_blockid_s1_idx].append( varpos ) 

                
            else:                                                  # A new block started, it can be the first block
                
                allele_blocks_s1.append( [int(allele_s1[0])] )     # alleles of a phase block (all vars are hetrozygous) 
                varpos_blocks_s1.append( [varpos] )                # variant posistion of a block
                blockid_blocks_s1.append( blockid_s1 )             # add the blockid of the first var  to the list of blockids
                
                previous_blockid_s1 = blockid_s1
                previous_blockid_s1_idx = len(blockid_blocks_s1)-1  # the index of the last block id
                


                
       
    
                                                ###  sample 2 ###
            split_s2 = s2.split(":")
            blockid_s2 = split_s2[idx_ps]
            allele_s2 = split_s2[idx_gt]
            if blockid_s2 == ".":
                pass
            elif blockid_s2 in blockid_blocks_s2:
                assert allele_s2[1] == "|", "allele_s2[1] should be | and so phased. A non-phased variant is in the phase block"
                
                if blockid_s2 != previous_blockid_s2:  
                    previous_blockid_s2_idx = blockid_blocks_s2.index(blockid_s2)
                    num_weird_occurs+=1
                allele_blocks_s2[previous_blockid_s2_idx].append( int(allele_s2[0]) )
                varpos_blocks_s2[previous_blockid_s2_idx].append( varpos ) 
                var_blockids1_blocks_s2[previous_blockid_s2_idx].append( blockid_s1 )
                
            else:                                                 # A new block started   (it can be the first block)
                allele_blocks_s2.append( [int(allele_s2[0])] )     # alleles of a phase block (all vars are hetrozygous) 
                varpos_blocks_s2.append( [varpos] )                # variant posistion of a block
                blockid_blocks_s2.append( blockid_s2 )             # add the blockid of the first var  to the list of blockids
                var_blockids1_blocks_s2.append( [blockid_s1] )  
                previous_blockid_s2 = blockid_s2
                previous_blockid_s2_idx = len(blockid_blocks_s2)-1   # the index of the last block id
                

                                                ###  sample 3 ###
            split_s3 = s3.split(":")
            blockid_s3 = split_s3[idx_ps]
            allele_s3 = split_s3[idx_gt]
            if blockid_s3 == ".":
                pass
            elif blockid_s3 in blockid_blocks_s3:
                assert allele_s3[1] == "|", "allele_s3[1] should be | and so phased. A non-phased variant is in the phase block"
                
                if blockid_s3 != previous_blockid_s3:  
                    previous_blockid_s3_idx = blockid_blocks_s3.index(blockid_s3)
                    num_weird_occurs+=1
                allele_blocks_s3[previous_blockid_s3_idx].append( int(allele_s3[0]) )
                varpos_blocks_s3[previous_blockid_s3_idx].append( varpos ) 
                var_blockids1_blocks_s3[previous_blockid_s3_idx].append( blockid_s1 )
                
            else:                                                 # A new block started   (it can be the first block)
                allele_blocks_s3.append( [int(allele_s3[0])] )     # alleles of a phase block (all vars are hetrozygous) 
                varpos_blocks_s3.append( [varpos] )                # variant posistion of a block
                blockid_blocks_s3.append( blockid_s3 )             # add the blockid of the first var  to the list of blockids
                var_blockids1_blocks_s3.append( [blockid_s1] )  
                previous_blockid_s3 = blockid_s3
                previous_blockid_s3_idx = len(blockid_blocks_s3)-1   # the index of the last block id
                
   
    
                                                ###  sample 3 ###
            split_s4 = s4.split(":")
            blockid_s4 = split_s4[idx_ps]
            allele_s4 = split_s4[idx_gt]
            if blockid_s4 == ".":
                pass
            elif blockid_s4 in blockid_blocks_s4:
                assert allele_s4[1] == "|", "allele_s4[1] should be | and so phased. A non-phased variant is in the phase block"
                
                if blockid_s4 != previous_blockid_s4:  
                    previous_blockid_s4_idx = blockid_blocks_s4.index(blockid_s4)
                    num_weird_occurs+=1
                allele_blocks_s4[previous_blockid_s4_idx].append( int(allele_s4[0]) )
                varpos_blocks_s4[previous_blockid_s4_idx].append( varpos ) 
                var_blockids1_blocks_s4[previous_blockid_s4_idx].append( blockid_s1 )
                
            else:                                                 # A new block started   (it can be the first block)
                allele_blocks_s4.append( [int(allele_s4[0])] )     # alleles of a phase block (all vars are hetrozygous) 
                varpos_blocks_s4.append( [varpos] )                # variant posistion of a block
                blockid_blocks_s4.append( blockid_s4 )             # add the blockid of the first var  to the list of blockids
                var_blockids1_blocks_s4.append( [blockid_s1] )  
                previous_blockid_s4 = blockid_s4
                previous_blockid_s4_idx = len(blockid_blocks_s4)-1   # the index of the last block id
                




print("Num phaseblocks in sample one is ",len(allele_blocks_s1))
print("Num phaseblocks in sample two is ",len(allele_blocks_s2))
print("Num phaseblocks in sample three is ",len(allele_blocks_s3))
print("Num phaseblocks in sample four is ",len(allele_blocks_s4))
print("Num vars in sample five-parental is ",len(varposs_s5))
print("Num all vars is ",len(varposs_all))

print("number of occurance of weird case (=two variants are in the same block, but far from each other) ", num_weird_occurs)
 


# In[2]:


# for test code
# cd /Volumes/work/myjupyter/phaseyou
#grep "46014176" 22_all.vcf

#blockid_blocks_s2.index("46014176")
#print(allele_blocks_s2[5568],"\n",varpos_blocks_s2[5568])


# In[ ]:





# In[47]:





# In[ ]:





# In[3]:


#%pwd
num_phased_s1=0
for varpos_block_s1 in varpos_blocks_s1:
    num_phased_s1+=len(varpos_block_s1)
    if len(varpos_block_s1)<2:
        print("small",varpos_block_s1)
print(num_phased_s1)


# num_phased_s2=0
# for allele_block_s2 in allele_blocks_s2:
#     num_phased_s2+=len(allele_block_s2)
#     if len(allele_block_s2)<2:
#         print("small,s2",len(allele_block_s2))
# print(num_phased_s2)

num_phased_a=0
for a in var_blockids1_blocks_s2:
    num_phased_a+=len(a)
print(num_phased_a)   
#len(var_blockids1_blocks_s2[100])
#print(var_blockids1_blocks_s2[0])


# In[ ]:





# In[ ]:





# In[4]:


# var_blockids1_blocks_s3
# for the first variant in the first block of s3 (varpos_blocks_s3[0][0] which is 16050159),  
# the same variant has the block id  of  var_blockids1_blocks_s3[0][0] (which is'16050159') in sample 1.


print(len(var_blockids1_blocks_s3),len(var_blockids1_blocks_s3[0]))
print(len(varpos_blocks_s3),len(varpos_blocks_s3[0]))


# In[5]:


intersected_blocks_idins1_s3=[] # for each block in s3, the corresponding s1 block id

for var_blockids1_blocks_s3_i in var_blockids1_blocks_s3:
    var_blockids1_blocks_s3_i_set=set(var_blockids1_blocks_s3_i)
    
    if '.' in var_blockids1_blocks_s3_i_set:
        var_blockids1_blocks_s3_i_set.discard('.')
        
    if len(var_blockids1_blocks_s3_i_set)==1: # all the variant in s3, corresponds to only one block in s1
        intersected_idins1_s3=list(var_blockids1_blocks_s3_i_set)[0]
        
    if len(var_blockids1_blocks_s3_i_set)==0: # the variants in s3 are not phased in s1
        intersected_idins1_s3=0
    if len(var_blockids1_blocks_s3_i_set)>1:
        # find the most frequent blockid in the phaseblock
        intersected_idins1_s3=max(var_blockids1_blocks_s3_i_set, key=var_blockids1_blocks_s3_i.count) 
    
        
    intersected_blocks_idins1_s3.append(intersected_idins1_s3)

    
    
    
print("number of blocks in s3 without any phased variant in s1 ",len([i for i in intersected_blocks_idins1_s3 if i==0]),"out of ", len(intersected_blocks_idins1_s3))
# now in the list  intersected_blocks_idins1_s3, for each block in s3, we know one block in s1. 


# In[ ]:





# In[ ]:





# In[ ]:





# In[6]:


number_phaseblocks_s3= len(allele_blocks_s3)
allele_blocks_s3_reordered = []

reordered_block_i_list_s3=[]
for block_i in range(number_phaseblocks_s3): # 
   
    allele_block_s3 = allele_blocks_s3[block_i]
    varpos_block_s3 = varpos_blocks_s3[block_i]

    intersected_idins1_s3 = intersected_blocks_idins1_s3[block_i]
    if intersected_idins1_s3 !=0:
        block_indx_s1 = blockid_blocks_s1.index(intersected_idins1_s3)
    varpos_block_s1 = varpos_blocks_s1[block_indx_s1]
    allele_block_s1 = allele_blocks_s1[block_indx_s1]
    
    intersect_varpos = set(varpos_block_s1) & set(varpos_block_s3)
    similiar_count = 0
    for varpos in intersect_varpos:
        allele_s1 = allele_block_s1[varpos_block_s1.index(varpos)]
        allele_s3 = allele_block_s3[varpos_block_s3.index(varpos)]
        if allele_s1 == allele_s3: similiar_count+=1
    #print(similiar_count/len(intersect_varpos))
    if similiar_count < len(intersect_varpos)/2:
        #print("need", similiar_count)
        allele_block_s3_reordered = [1-i for i in allele_block_s3]
        #allele_blocks_s3[block_i] = allele_block_s3_reordered
        allele_blocks_s3_reordered.append(allele_block_s3_reordered)
        #print("it needs to re-order")
        reordered_block_i_list_s3.append(block_i)
    else:

        allele_blocks_s3_reordered.append(allele_block_s3)

print("numebr of re-ordered block s3 in accordance with s1", len(reordered_block_i_list_s3),"out of ",len(allele_blocks_s3_reordered))


# In[ ]:





# In[7]:


# for test
#for i in [0,1,2,2000]:
#    print(i in reordered_block_i_list,allele_blocks_s3[i],allele_blocks_s3_reordered[i],varpos_blocks_s3[i])
## grep -v "#" 22_5sample.vcf | head -n 10 | cut -f 2,10,12


# In[8]:


intersected_blocks_idins1_s4=[] # for each block in s4, the corresponding s1 block id

for var_blockids1_blocks_s4_i in var_blockids1_blocks_s4:
    var_blockids1_blocks_s4_i_set=set(var_blockids1_blocks_s4_i)
    
    if '.' in var_blockids1_blocks_s4_i_set:
        var_blockids1_blocks_s4_i_set.discard('.')
        
    if len(var_blockids1_blocks_s4_i_set)==1:
        intersected_idins1_s4=list(var_blockids1_blocks_s4_i_set)[0]
    if len(var_blockids1_blocks_s4_i_set)==0:
        intersected_idins1_s4=0
    if len(var_blockids1_blocks_s4_i_set)>1:
        # find the most frequent blockid in the phaseblock
        intersected_idins1_s4=max(var_blockids1_blocks_s4_i_set, key=var_blockids1_blocks_s4_i.count) 
    
        
    intersected_blocks_idins1_s4.append(intersected_idins1_s4)


# In[9]:


number_phaseblocks_s4= len(allele_blocks_s4)
allele_blocks_s4_reordered = []

reordered_block_i_list_s4=[]
for block_i in range(number_phaseblocks_s4): # 
   
    allele_block_s4 = allele_blocks_s4[block_i]
    varpos_block_s4 = varpos_blocks_s4[block_i]

    intersected_idins1_s4 = intersected_blocks_idins1_s4[block_i]
    if intersected_idins1_s4 !=0:
        block_indx_s1 = blockid_blocks_s1.index(intersected_idins1_s4)
    varpos_block_s1 = varpos_blocks_s1[block_indx_s1]
    allele_block_s1 = allele_blocks_s1[block_indx_s1]
    
    intersect_varpos = set(varpos_block_s1) & set(varpos_block_s4)
    similiar_count = 0
    for varpos in intersect_varpos:
        allele_s1 = allele_block_s1[varpos_block_s1.index(varpos)]
        allele_s4 = allele_block_s4[varpos_block_s4.index(varpos)]
        if allele_s1 == allele_s4: similiar_count+=1
    #print(similiar_count/len(intersect_varpos))
    if similiar_count < len(intersect_varpos)/2:
        #print("need", similiar_count)
        allele_block_s4_reordered = [1-i for i in allele_block_s4]
        #allele_blocks_s4[block_i] = allele_block_s4_reordered
        allele_blocks_s4_reordered.append(allele_block_s4_reordered)
        #print("it needs to re-order")
        reordered_block_i_list_s4.append(block_i)
    else:

        allele_blocks_s4_reordered.append(allele_block_s4)

print("numebr of re-ordered block s4 in accordance with s1", len(reordered_block_i_list_s4),"out of ",len(allele_blocks_s4_reordered))


# In[10]:


## for test
# for i in [0,1,2,2000]:
#     print(i in reordered_block_i_list,allele_blocks_s4[i],allele_blocks_s4_reordered[i],varpos_blocks_s4[i])

## grep -v "#" 22_5sample.vcf | head -n 10 | cut -f 2,10,12
## grep 24247293 22_5sample.vcf  | head -n 10 | cut -f 2,10,13 


# In[ ]:





# In[11]:


intersected_blocks_idins1_s2=[] # for each block in s2, the corresponding block_id of s1
ratio_all = []
for var_blockids1_blocks_s2_i in var_blockids1_blocks_s2:
    var_blockids1_blocks_s2_i_set=set(var_blockids1_blocks_s2_i)
    
    if '.' in var_blockids1_blocks_s2_i_set:
        var_blockids1_blocks_s2_i_set.discard('.')
        
    if len(var_blockids1_blocks_s2_i_set)==1:
        intersected_idins1_s2=list(var_blockids1_blocks_s2_i_set)[0]
        ratio = 1
        
    if len(var_blockids1_blocks_s2_i_set)==0:
        intersected_idins1_s2=0
    if len(var_blockids1_blocks_s2_i_set)>1:
        # find the most frequent blockid in the phaseblock
        intersected_idins1_s2=max(var_blockids1_blocks_s2_i_set, key=var_blockids1_blocks_s2_i.count) 

        ratio=var_blockids1_blocks_s2_i.count(intersected_idins1_s2)/len(var_blockids1_blocks_s2_i)
        if ratio<.71: print(var_blockids1_blocks_s2_i[0])
    ratio_all.append(ratio)
    
        
    intersected_blocks_idins1_s2.append(intersected_idins1_s2)


# In[ ]:





# In[ ]:





# In[12]:


print(sum(ratio_all)/len(ratio_all))
print(min(ratio_all))


# In[13]:


import numpy as np
np.sort([round(i,2) for i in ratio_all if i < .97])


# In[52]:


number_phaseblocks_s2= len(allele_blocks_s2)
allele_blocks_s2_reordered = []

reordered_block_i_list_s2=[]
for block_i in range(number_phaseblocks_s2): # 
   
    allele_block_s2 = allele_blocks_s2[block_i]
    varpos_block_s2 = varpos_blocks_s2[block_i]

    intersected_idins1_s2 = intersected_blocks_idins1_s2[block_i]
    if intersected_idins1_s2 !=0:
        block_indx_s1 = blockid_blocks_s1.index(intersected_idins1_s2)
    varpos_block_s1 = varpos_blocks_s1[block_indx_s1]
    allele_block_s1 = allele_blocks_s1[block_indx_s1]
    
    intersect_varpos = set(varpos_block_s1) & set(varpos_block_s2)
    similiar_count = 0
    for varpos in intersect_varpos:
        allele_s1 = allele_block_s1[varpos_block_s1.index(varpos)]
        allele_s2 = allele_block_s2[varpos_block_s2.index(varpos)]
        if allele_s1 == allele_s2: similiar_count+=1
    #print(similiar_count/len(intersect_varpos))
    if block_i==17:
        print(similiar_count,len(intersect_varpos),len(varpos_block_s2),len(varpos_block_s1))
    
    if similiar_count < len(intersect_varpos)/2:
        #print("need", similiar_count)
        allele_block_s2_reordered = [1-i for i in allele_block_s2]
        #allele_blocks_s2[block_i] = allele_block_s2_reordered
        allele_blocks_s2_reordered.append(allele_block_s2_reordered)
        #print("it needs to re-order")
        reordered_block_i_list_s2.append(block_i)
    else:

        allele_blocks_s2_reordered.append(allele_block_s2)

print("numebr of re-ordered block s2 in accordance with s1", len(reordered_block_i_list_s2),"out of ",len(allele_blocks_s2_reordered))


# In[51]:


1249 <2340/2


# In[ ]:





# In[15]:


matrix_row_lists=[]
for varpos in varposs_all:
    # for increasing the speed we can use the information that the phase blocks are in order.
    # but not for those weird phase block


    allele_s1='.'
    blockid_s1='.'
    for block_i_s1 in range(len(blockid_blocks_s1)):
        allele_block_s1 = allele_blocks_s1[block_i_s1]
        varpos_block_s1 = varpos_blocks_s1[block_i_s1]

        if varpos in varpos_block_s1:
            allele_s1 = allele_block_s1[varpos_block_s1.index(varpos)]
            blockid_s1 = blockid_blocks_s1[block_i_s1]
            break # when var is found, no need to continue for loop

            
    allele_s2='.'
    blockid_s2='.'
    for block_i_s2 in range(len(blockid_blocks_s2)):
        allele_block_s2 = allele_blocks_s2_reordered[block_i_s2]
        varpos_block_s2 = varpos_blocks_s2[block_i_s2]

        if varpos in varpos_block_s2:
            allele_s2 = allele_block_s2[varpos_block_s2.index(varpos)]
            blockid_s2 = blockid_blocks_s2[block_i_s2]
            break # when var is found, no need to continue for loop

    allele_s3='.'
    blockid_s3='.'
    for block_i_s3 in range(len(blockid_blocks_s3)):
        allele_block_s3 = allele_blocks_s3_reordered[block_i_s3]
        varpos_block_s3 = varpos_blocks_s3[block_i_s3]

        if varpos in varpos_block_s3:
            allele_s3 = allele_block_s3[varpos_block_s3.index(varpos)]
            blockid_s3 = blockid_blocks_s3[block_i_s3]
            break # when var is found, no need to continue for loop

    allele_s4='.'
    blockid_s4='.'
    for block_i_s4 in range(len(blockid_blocks_s4)):
        allele_block_s4 = allele_blocks_s4_reordered[block_i_s4]
        varpos_block_s4 = varpos_blocks_s4[block_i_s4]

        if varpos in varpos_block_s4:
            allele_s4 = allele_block_s4[varpos_block_s4.index(varpos)]
            blockid_s4 = blockid_blocks_s4[block_i_s4]
            break # when var is found, no need to continue for loop

    
            
    if varpos in varposs_s5:
        allele_s5 = alleles_s5[varposs_s5.index(varpos)]
    else: allele_s5= '.'

    if not (allele_s1=='.' and allele_s2=='.' and allele_s3=='.' and allele_s4=='.' and allele_s5=='.'): 
        matrix_row_list=[str(varpos),str(allele_s1)+":"+str(blockid_s1),
                         str(allele_s2)+":"+str(blockid_s2),str(allele_s3)+":"+str(blockid_s3),
                         str(allele_s4)+":"+str(blockid_s4), str(allele_s5)] 
        matrix_row_lists.append(matrix_row_list)

        
# one option can be for loop over blockid_blocks_s1 instead of varposs_all, then search in  blockid_blocks_s1 also for printing sort them to overcome weird phase block
# the issue will be for those variant that are not phased in blockid_blocks_s1 but in blockid_blocks_s2


# In[ ]:





# In[ ]:





# In[19]:


matrix_row_lists[0]


# In[20]:


matrix_row_lists[0]



# In[21]:



import copy 


# In[23]:


matrix_row_lists_maj = copy.deepcopy(matrix_row_lists)

#matrix_row_lists_maj = copy.deepcopy(matrix_row_lists) # matrix_row_lists[:]
num_err_maj=0

num_shared_majv=0
num_shared_s1=0
num_shared_s2=0
num_shared_s3=0
num_shared_s4=0

num_err_s1=0
num_err_s2=0
num_err_s3=0
num_err_s4=0
num_maj=0
num1,num2,num3,num4=[0,0,0,0]




for row_i, matrix_row_list in enumerate(matrix_row_lists):
    
    varpos, s1, s2, s3, s4, allele_s5 = matrix_row_list

    
    all_alleles_raw=[s1[0], s2[0],s3[0],s4[0]]
    all_alleles=[i for i in all_alleles_raw if i != '.']   
    if all_alleles:
        if all_alleles.count('0') > all_alleles.count('1'):
            allele_majvt = '0'
        elif all_alleles.count('1') > all_alleles.count('0'):
            allele_majvt = '1'
        elif all_alleles.count('1') == all_alleles.count('0'):
            allele_majvt = all_alleles[0] # s1

        # allele_majvt=max(set(all_alleles), key=all_alleles.count)
    else:
        allele_majvt="."

    if s1==".":
        blockid_s1="."
    else:
        blockid_s1= s1.split(":")[1]

    allele_majvt=allele_majvt+":"+blockid_s1

    
    
    
    matrix_row_lists_maj[row_i].append(allele_majvt)
    
    
    if allele_majvt!='.':
        num_maj+=1
        if allele_s5 !='.':
            num_shared_majv+=1
            if allele_s5!=allele_majvt:
                num_err_maj+=1
    
    if s1[0] !='.':
        num1+=1
        if allele_s5 !='.' :
            num_shared_s1 +=1
            if allele_s5!=s1[0]:
                num_err_s1+=1

    
    if s2[0] !='.':
        num2+=1
        if allele_s5 !='.' :
            num_shared_s2 +=1
            if allele_s5!=s2[0]:
                num_err_s2+=1
                
    
    if s3[0] !='.':
        num3+=1
        if allele_s5 !='.' :
            num_shared_s3 +=1
            if allele_s5!=s3[0]:
                num_err_s3+=1
                
    
    if s4[0] !='.':
        num4+=1
        if allele_s5 !='.' :
            num_shared_s4 +=1
            if allele_s5!=s4[0]:
                num_err_s4+=1
            
print("Majoritvote\t",num_maj , num_shared_majv, num_err_maj, round(100*num_err_maj/num_shared_majv,2))
print("hifi-whatshp \t",num1, num_shared_s1,num_err_s1, round(100*num_err_s1/num_shared_s1,2))
print("hifi-hapcut \t",num2, num_shared_s2,num_err_s2, round(100*num_err_s2/num_shared_s2,2))
print("illum-whatshp  \t", num3, num_shared_s3,num_err_s3, round(100*num_err_s3/num_shared_s3,2))
print("illum-hapcut \t",num4, num_shared_s4,num_err_s4, round(100*num_err_s4/num_shared_s4,2))

  
# ? ?? s1 is not reordered based on parental ???


# In[24]:


max(set(all_alleles), key=all_alleles.count)
all_alleles

all_alleles.count('0')


# In[ ]:





# In[26]:


matrix_row_list=matrix_row_lists[10]
matrix_row_list
varpos, s1, s2, s3, s4, allele_s5 = matrix_row_list


# In[ ]:




if s1==".":
    blockid_s1="."
else:
    blockid_s1= s1.split(":")[1]
                            
#blockid_s1

matrix_row_list_prd


# In[ ]:


if s1==".":
    blockid_s1="."
else:
    blockid_s1= s1.split(":")[1]
matrix_row_list_prd[1].split(":")[1]


# In[ ]:





# In[27]:


diff_list=[]

matrix_row_lists_prd= var_lines_list
var_lines_list_prd=[]

varpos_list_inmatrix=[i[0] for i in matrix_row_lists_maj]

for var_line_list in var_lines_list: # 
    var_line_list_part=var_line_list[:9]

    varpos = var_line_list_part[1]
    if varpos in varpos_list_inmatrix:
        
        gt_flags = var_line_list[8]
        idx_ps = gt_flags.split(":").index("PS")
        idx_gt = gt_flags.split(":").index("GT")   
        
        var_idx = varpos_list_inmatrix.index(varpos)        
        matrix_row_list_prd = matrix_row_lists_maj[var_idx]  
        
        allele_prd_blockId_prd = matrix_row_list_prd[6]
        #allele_prd_blockId_prd = matrix_row_list_prd[7]
        allele_prd, blockId_prd=allele_prd_blockId_prd.split(":")
        blockId_s1 = matrix_row_list_prd[1].split(":")[1]
        
        
        if allele_prd != "." and blockId_prd != ".":
            alleles_prd = allele_prd+"|"+str(1-int(allele_prd))+":"+blockId_prd
        else:
            alleles_prd = "./.:."
        
#         if allele_prd != "." and blockId_prd == "." and blockId_s1_previous!= ".":
#             alleles_prd = allele_prd+"|"+str(1-int(allele_prd))+":"+blockId_s1       
        
        if blockId_prd != "." and  allele_prd != matrix_row_list_prd[1].split(":")[0]:
            diff_list.append((varpos,allele_prd,matrix_row_list_prd[1].split(":")[0]))
    
        
        var_line_list_prd=var_line_list[:8] +['GT:PS', alleles_prd]        
        var_lines_list_prd.append(var_line_list_prd)
        
        blockId_s1_previous =  blockId_s1
        
        


   


# In[28]:


len(diff_list)


# In[29]:


diff_list[:10]


# In[48]:


correct =0
wrong = 0
other = 0
for i in diff_list:
    var_pos=i[0]
    var_idx = varpos_list_inmatrix.index(var_pos)        
    matrix_row_list_prd = matrix_row_lists_maj[var_idx] 
    if  matrix_row_list_prd[5]!='.':
        
        if matrix_row_list_prd[5]==matrix_row_list_prd[6][0] and matrix_row_list_prd[1][0]!=matrix_row_list_prd[6][0]:
            correct+=1
            print("correct_prediction",matrix_row_list_prd)
        elif matrix_row_list_prd[5]!=matrix_row_list_prd[6][0] and matrix_row_list_prd[1][0]!=matrix_row_list_prd[6][0]:
            wrong+=1
            #print("wrong_prediction",matrix_row_list_prd)
    else:
        other+=1
            
print(correct,wrong,other)


# In[ ]:





# In[ ]:





# In[31]:



vcf_out_address=vcf_file_address+"_majority_3b.vcf"

vcf_out_file=open(vcf_out_address,'w');

for header_line in header_lines_list[:-1]: # except last line should be edit for one sample not five
    vcf_out_file.write(header_line+"\n")

vcf_out_file.write(header_lines_list[-1][:51]+"\n")
    
for var_line_list in var_lines_list_prd:
    
    vcf_out_file.write("\t".join(var_line_list)+"\n")

vcf_out_file.close()


# In[ ]:


#print(blockid_blocks_s1)
#for i in range(len(allele_blocks_s1)):
#    print(len(allele_blocks_s1[i]))
    
#print(blockid_blocks_s2)
#for i in range(len(allele_blocks_s2)):
#    print(len(allele_blocks_s2[i]))


# In[ ]:





# In[ ]:



# # re-ordering  the phase block allele based on parent

# ## important ## 
# # if we don't have the parental , how we can reorder


# number_phaseblocks_s1= len(allele_blocks_s1)

# #allele_blocks_s1_new=list(allele_blocks_s1)
# for block_i in range(number_phaseblocks_s1): # 
   
#     allele_block_s1=allele_blocks_s1[block_i]
#     varpos_block_s1=varpos_blocks_s1[block_i]


#     similiar_count=0
#     diss_count=0
#     for var_i, varpos in enumerate(varpos_block_s1):
#         allele_s1= allele_block_s1[var_i]
        
#         if varpos in varposs_s5:
#             allele_s5= alleles_s5[varposs_s5.index(varpos)]
#             if allele_s1 == allele_s5:
#                 similiar_count+=1
#             else: 
#                 diss_count+=1
                
#     #print(similiar_count, diss_count)          
#     if similiar_count < diss_count:
#         allele_block_s1_reordered = [1-i for i in allele_block_s1] # converting 0->1 and 1->0
#         allele_blocks_s1[block_i]=allele_block_s1_reordered
#         #print("it needs to re-order")
    
# #         allele_block_s1_new = [1-i for i in allele_block_s1]
# #         allele_blocks_s1_new[block_i]=allele_block_s1_new
        


# In[ ]:



# # re-ordering  the phase block allel based on parent

# ## important ?? 
# # if we don't have the parental , how we can reorder


# number_phaseblocks_s2= len(allele_blocks_s2)
# #allele_blocks_s1_new=list(allele_blocks_s1)
# for block_i in range(number_phaseblocks_s2): # 
   
#     allele_block_s2=allele_blocks_s2[block_i]
#     varpos_block_s2=varpos_blocks_s2[block_i]


#     similiar_count=0
#     diss_count=0
#     for var_i, varpos in enumerate(varpos_block_s2):
#         allele_s2= allele_block_s2[var_i]
        
#         if varpos in varposs_s5:
#             allele_s5= alleles_s5[varposs_s5.index(varpos)]
#             if allele_s2 == allele_s5:
#                 similiar_count+=1
#             else: 
#                 diss_count+=1
                
#     #print(similiar_count, diss_count)          
#     if similiar_count < diss_count:
#         allele_block_s2_reordered = [1-i for i in allele_block_s2] # converting 0->1 and 1->0
#         allele_blocks_s2[block_i]=allele_block_s2_reordered


        
# print("s2 done")     
# number_phaseblocks_s3= len(allele_blocks_s3)
# for block_i in range(number_phaseblocks_s3): # 
   
#     allele_block_s3=allele_blocks_s3[block_i]
#     varpos_block_s3=varpos_blocks_s3[block_i]


#     similiar_count=0
#     diss_count=0
#     for var_i, varpos in enumerate(varpos_block_s3):
#         allele_s3= allele_block_s3[var_i]
        
#         if varpos in varposs_s5:
#             allele_s5= alleles_s5[varposs_s5.index(varpos)]
#             if allele_s3 == allele_s5:
#                 similiar_count+=1
#             else: 
#                 diss_count+=1
#     #print(similiar_count, diss_count)          
#     if similiar_count < diss_count:
#         allele_block_s3_reordered = [1-i for i in allele_block_s3] # converting 0->1 and 1->0
#         allele_blocks_s3[block_i]=allele_block_s3_reordered
    

# print("s3 done")     
# number_phaseblocks_s4= len(allele_blocks_s4)
# for block_i in range(number_phaseblocks_s4): # 
   
#     allele_block_s4=allele_blocks_s4[block_i]
#     varpos_block_s4=varpos_blocks_s4[block_i]


#     similiar_count=0
#     diss_count=0
#     for var_i, varpos in enumerate(varpos_block_s4):
#         allele_s4= allele_block_s4[var_i]
        
#         if varpos in varposs_s5:
#             allele_s5= alleles_s5[varposs_s5.index(varpos)]
#             if allele_s4 == allele_s5:
#                 similiar_count+=1
#             else: 
#                 diss_count+=1
                
#     #print(similiar_count, diss_count)          
#     if similiar_count < diss_count:
#         allele_block_s4_reordered = [1-i for i in allele_block_s4] # converting 0->1 and 1->0
#         allele_blocks_s4[block_i]=allele_block_s4_reordered
    


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


# number_phaseblocks_s2= len(allele_blocks_s2)

# for block_i in range(number_phaseblocks_s2): # 
   
#     allele_block_s2=allele_blocks_s2[block_i]
#     varpos_block_s2=varpos_blocks_s2[block_i]

#     intersected_idins1_s2=intersected_blocks_idins1_s2[block_i]
#     if intersected_idins1_s2 !=0:
#         block_indx_s1= blockid_blocks_s1.index(intersected_idins1_s2)
#     varpos_block_s1=varpos_blocks_s1[block_indx_s1]
#     allele_block_s1=allele_blocks_s1[block_indx_s1]
    
#     intersect_varpos=set(varpos_block_s1) & set(varpos_block_s2)
#     similiar_count=0
#     for varpos in intersect_varpos:
#         allele_s1= allele_block_s1[varpos_block_s1.index(varpos)]
#         allele_s2= allele_block_s2[varpos_block_s2.index(varpos)]
#         if allele_s1==allele_s2: similiar_count+=1
#     #print(similiar_count/len(intersect_varpos))
#     if similiar_count < len(intersect_varpos)/2:
#         #print("need", similiar_count)
#         allele_block_s2_reordered= [1-i for i in allele_block_s2]
#         allele_blocks_s2[block_i]=allele_block_s2_reordered
#         #print("it needs to re-order")


# In[ ]:



        
# number_phaseblocks_s4= len(allele_blocks_s4)

# for block_i in range(number_phaseblocks_s4): # 
   
#     allele_block_s4=allele_blocks_s4[block_i]
#     varpos_block_s4=varpos_blocks_s4[block_i]

#     intersected_idins1_s4=intersected_blocks_idins1_s4[block_i]
#     if intersected_idins1_s4 !=0:
#         block_indx_s1= blockid_blocks_s1.index(intersected_idins1_s4)
#     varpos_block_s1=varpos_blocks_s1[block_indx_s1]
#     allele_block_s1=allele_blocks_s1[block_indx_s1]
    
#     intersect_varpos=set(varpos_block_s1) & set(varpos_block_s4)
#     similiar_count=0
#     for varpos in intersect_varpos:
#         allele_s1= allele_block_s1[varpos_block_s1.index(varpos)]
#         allele_s4= allele_block_s4[varpos_block_s4.index(varpos)]
#         if allele_s1==allele_s4: similiar_count+=1
#     #print(similiar_count/len(intersect_varpos))
#     if similiar_count < len(intersect_varpos)/2:
#         #print("need", similiar_count)
#         allele_block_s4_reordered= [1-i for i in allele_block_s4]
#         allele_blocks_s4[block_i]=allele_block_s4_reordered
#         #print("it needs to re-order")


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


# matrix_row_lists=[]
# for varpos in varposs_all:
#     # for increasing the speed we can use the information that the phase blocks are in order.
#     # but not for those weird phase block


#     allele_s1='.'
#     blockid_s1='.'
#     for block_i_s1 in range(len(blockid_blocks_s1)):
#         allele_block_s1 = allele_blocks_s1[block_i_s1]
#         varpos_block_s1 = varpos_blocks_s1[block_i_s1]

#         if varpos in varpos_block_s1:
#             allele_s1 = allele_block_s1[varpos_block_s1.index(varpos)]
#             blockid_s1 = blockid_blocks_s1[block_i_s1]
#             break # when var is found, no need to continue for loop

            
#     allele_s2='.'
#     blockid_s2='.'
#     for block_i_s2 in range(len(blockid_blocks_s2)):
#         allele_block_s2 = allele_blocks_s2[block_i_s2]
#         varpos_block_s2 = varpos_blocks_s2[block_i_s2]

#         if varpos in varpos_block_s2:
#             allele_s2 = allele_block_s2[varpos_block_s2.index(varpos)]
#             blockid_s2 = blockid_blocks_s2[block_i_s2]
#             break # when var is found, no need to continue for loop

#     allele_s3='.'
#     blockid_s3='.'
#     for block_i_s3 in range(len(blockid_blocks_s3)):
#         allele_block_s3 = allele_blocks_s3[block_i_s3]
#         varpos_block_s3 = varpos_blocks_s3[block_i_s3]

#         if varpos in varpos_block_s3:
#             allele_s3 = allele_block_s3[varpos_block_s3.index(varpos)]
#             blockid_s3 = blockid_blocks_s3[block_i_s3]
#             break # when var is found, no need to continue for loop

#     allele_s4='.'
#     blockid_s4='.'
#     for block_i_s4 in range(len(blockid_blocks_s4)):
#         allele_block_s4 = allele_blocks_s4[block_i_s4]
#         varpos_block_s4 = varpos_blocks_s4[block_i_s4]

#         if varpos in varpos_block_s4:
#             allele_s4 = allele_block_s4[varpos_block_s4.index(varpos)]
#             blockid_s4 = blockid_blocks_s4[block_i_s4]
#             break # when var is found, no need to continue for loop

    
            
#     if varpos in varposs_s5:
#         allele_s5 = alleles_s5[varposs_s5.index(varpos)]
#     else: allele_s5= '.'

#     if not (allele_s1=='.' and allele_s2=='.' and allele_s3=='.' and allele_s4=='.' and allele_s5=='.'): 
#         matrix_row_list=[str(varpos),str(allele_s1)+":"+str(blockid_s1),
#                          str(allele_s2)+":"+str(blockid_s2),str(allele_s3)+":"+str(blockid_s3),
#                          str(allele_s4)+":"+str(blockid_s4), str(allele_s5)] 
#         matrix_row_lists.append(matrix_row_list)

        
# # one option can be for loop over blockid_blocks_s1 instead of varposs_all, then search in  blockid_blocks_s1 also for printing sort them to overcome weird phase block
# # the issue will be for those variant that are not phased in blockid_blocks_s1 but in blockid_blocks_s2


# In[ ]:


matrix_row_list


# In[ ]:


# matrix_row_lists_maj =matrix_row_lists[:]
# num_err_maj=0

# num_shared_majv=0
# num_shared_s1=0
# num_shared_s2=0
# num_shared_s3=0
# num_shared_s4=0

# num_err_s1=0
# num_err_s2=0
# num_err_s3=0
# num_err_s4=0
# num_maj=0
# num1,num2,num3,num4=[0,0,0,0]




# for row_i, matrix_row_list in enumerate(matrix_row_lists):
    
#     varpos, s1, s2, s3, s4, allele_s5 = matrix_row_list

    
#     all_alleles_raw=[s1[0], s2[0],s3[0],s4[0]]
#     all_alleles=[i for i in all_alleles_raw if i != '.']   
#     if all_alleles:
#         allele_majvt=max(set(all_alleles), key=all_alleles.count)
#     else:
#         allele_majvt="."
        
#     matrix_row_lists_maj[row_i].append(allele_majvt)
    
    
#     if allele_majvt!='.':
#         num_maj+=1
#         if allele_s5 !='.':
#             num_shared_majv+=1
#             if allele_s5!=allele_majvt:
#                 num_err_maj+=1
    
#     if s1[0] !='.':
#         num1+=1
#         if allele_s5 !='.' :
#             num_shared_s1 +=1
#             if allele_s5!=s1[0]:
#                 num_err_s1+=1

    
#     if s2[0] !='.':
#         num2+=1
#         if allele_s5 !='.' :
#             num_shared_s2 +=1
#             if allele_s5!=s2[0]:
#                 num_err_s2+=1
                
    
#     if s3[0] !='.':
#         num3+=1
#         if allele_s5 !='.' :
#             num_shared_s3 +=1
#             if allele_s5!=s3[0]:
#                 num_err_s3+=1
                
    
#     if s4[0] !='.':
#         num4+=1
#         if allele_s5 !='.' :
#             num_shared_s4 +=1
#             if allele_s5!=s4[0]:
#                 num_err_s4+=1
            
# print("Majoritvote\t",num_maj , num_shared_majv, num_err_maj, round(100*num_err_maj/num_shared_majv,2))
# print("hifi-whatshp \t",num1, num_shared_s1,num_err_s1, round(100*num_err_s1/num_shared_s1,2))
# print("hifi-hapcut \t",num2, num_shared_s2,num_err_s2, round(100*num_err_s2/num_shared_s2,2))
# print("illum-whatshp  \t", num3, num_shared_s3,num_err_s3, round(100*num_err_s3/num_shared_s3,2))
# print("illum-hapcut \t",num4, num_shared_s4,num_err_s4, round(100*num_err_s4/num_shared_s4,2))

  


# In[ ]:


print(matrix_row_lists[0])
var_line_list_part[1]


# In[ ]:


var_idx = varpos_list_inmatrix.index(varpos)        
matrix_row_list_prd = matrix_row_lists_maj[var_idx]   
matrix_row_list_prd


allele_prd = matrix_row_list_prd[6]


var_line_list_part.append(allele_prd+"|"+str(1-int(allele_prd))+":1")
#var_line_list_part



# In[ ]:


matrix_row_lists_prd= var_lines_list
var_lines_list_prd=[]

varpos_list_inmatrix=[i[0] for i in matrix_row_lists_maj]

for var_line_list in var_lines_list: # 
    var_line_list_part=var_line_list[:9]

    varpos = var_line_list_part[1]
    if varpos in varpos_list_inmatrix:
        
        gt_flags = var_line_list[8]
        idx_ps = gt_flags.split(":").index("PS")
        idx_gt = gt_flags.split(":").index("GT")   
        
        var_idx = varpos_list_inmatrix.index(varpos)        
        matrix_row_list_prd = matrix_row_lists_maj[var_idx]  
        

        allele_prd = matrix_row_list_prd[6]
        
        if allele_prd != ".":
            alleles_prd = allele_prd+"|"+str(1-int(allele_prd))+":1"
        else:
            alleles_prd = "./.:."
        
        
        
        var_line_list_prd=var_line_list[:8] +['GT:PS', alleles_prd]        
        var_lines_list_prd.append(var_line_list_prd)
        
        



vcf_out_address=vcf_file_address+"_majority.vcf"

vcf_out_file=open(vcf_out_address,'w');

for header_line in header_lines_list: # except last line should be edit for one sample not five
    vcf_out_file.write(header_line+"\n")
    
for var_line_list in var_lines_list_prd:
    
    vcf_out_file.write("\t".join(var_line_list)+"\n")

vcf_out_file.close()
   


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


matrix_address=vcf_file_address+"_matrix5_prd3.txt"

matrix_file=open(matrix_address,'w');
matrix_file.write("\t\t".join(["Pos","HiFiWsp","HiFiHpct","IllmWsp","IllmHpct","parental","predicted"])+"\n") # "Pos","ONT","CLR","HiFi","Illm","parental"
for matrix_row_list in matrix_row_lists_maj:

        matrix_file.write("\t".join(matrix_row_list)+"\n")

matrix_file.close()
   


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


result_consensus =[]

changes_varpos_pre=[]
changes_varpos_next=[]

matrix_row_lists_prd =[]
for row_i, matrix_row_list in enumerate(matrix_row_lists):
    varpos, s1, s2, s3, s4, allele_s5 = matrix_row_list
    
    allele_s1='.';blockid_s1='.'
    allele_s2='.';blockid_s2='.'
    allele_s3='.';blockid_s3='.'
    allele_s4='.';blockid_s4='.'
    if s1 != '.:.':
        allele_s1 = s1[0]
        blockid_s1 = s1[2:]
        blockid_s1_last_known = blockid_s1
        
    if s2 != '.:.':
        allele_s2 = s2[0]
        blockid_s2 = s2[2:]        
    if s3 != '.:.':
        allele_s3 = s3[0]
        blockid_s3 = s3[2:]
    if s4 != '.:.':
        allele_s4 = s4[0]
        blockid_s4 = s4[2:]
    
    allele_s_p = allele_s1
    blockid_s_p = blockid_s1
    
    if allele_s1 == '.' and allele_s2 != ".":

        if blockid_s2 == blockid_s1_last_known:
            allele_s_p = allele_s2
            blockid_s_p = blockid_s2
            changes_varpos_pre.append(varpos)
        else:             # check next variant, 

            blockid_s1_next = matrix_row_lists[row_i+1][1][2:]
            blockid_s2_next = matrix_row_lists[row_i+1][2][2:]
            if blockid_s2_next == blockid_s2:
                allele_s_p = allele_s2
                blockid_s_p = blockid_s2 
                changes_varpos_next.append(varpos)

            else:
                print("blockIds of s1 and s2 don't match: ",varpos,allele_s_p,blockid_s_p,)
            
    if allele_s1 != '.' :
        if allele_s1 != allele_s2 and allele_s2 == allele_s3 and allele_s3 == allele_s4:
            result_consensus.append(allele_s1 == allele_s5)
            #print(varpos, allele_s1 == allele_s5)
            allele_s_p = str(1-int(allele_s_p))

   
            
            
        
    matrix_row_list_prd = [str(varpos),str(allele_s1)+":"+str(blockid_s1),
                                  str(allele_s2)+":"+str(blockid_s2),str(allele_s3)+":"+str(blockid_s3),
                                  str(allele_s4)+":"+str(blockid_s4), str(allele_s5),
                                  str(allele_s_p)+":"+str(blockid_s_p)] 
    
    matrix_row_lists_prd.append(matrix_row_list_prd)
    
    
    


# In[ ]:


#sum(result_consensus), len(result_consensus)


# In[ ]:





# In[ ]:


#print(len(changes_varpos_pre),len(changes_varpos_next))


# In[ ]:


#changes_varpos_next[:20]


# In[ ]:





# In[ ]:


matrix_address=vcf_file_address+"_matrix5_p2.txt"

matrix_file=open(matrix_address,'w');
matrix_file.write("\t\t".join(["Pos","HiFiWsp","HiFiHpct","IllmWsp","IllmHpct","parental","predicted"])+"\n") # "Pos","ONT","CLR","HiFi","Illm","parental"
for matrix_row_list in matrix_row_lists_prd:

        matrix_file.write("\t".join(matrix_row_list)+"\n")

matrix_file.close()
   


# In[ ]:


vcf_file_address


# In[ ]:



result_consensus =[]

changes_varpos_pre=[]
changes_varpos_next=[]

matrix_row_lists_prd =[]


cntr1=cntr2=cntr3=cntr4=cntr0=0
for row_i, matrix_row_list in enumerate(matrix_row_lists):
    varpos, s1, s2, s3, s4, allele_s5 = matrix_row_list
    
    allele_s1='.';blockid_s1='.'
    allele_s2='.';blockid_s2='.'
    allele_s3='.';blockid_s3='.'
    allele_s4='.';blockid_s4='.'
    if s1 != '.:.':
        allele_s1 = s1[0]
        blockid_s1 = s1[2:]
        blockid_s1_last_known = blockid_s1
        
    if s2 != '.:.':
        allele_s2 = s2[0]
        blockid_s2 = s2[2:]        
    if s3 != '.:.':
        allele_s3 = s3[0]
        blockid_s3 = s3[2:]
    if s4 != '.:.':
        allele_s4 = s4[0]
        blockid_s4 = s4[2:]
    
    allele_s_p = allele_s1
    blockid_s_p = blockid_s1
    
    if allele_s1 == '.' and allele_s2 != ".":
        
        

        if blockid_s2 == blockid_s1_last_known:
            allele_s_p = allele_s2
            blockid_s_p = blockid_s2
            changes_varpos_pre.append(varpos)
        else:             # check next variant, 

            blockid_s1_next = matrix_row_lists[row_i+1][1][2:]
            blockid_s2_next = matrix_row_lists[row_i+1][2][2:]
            if blockid_s2_next == blockid_s2:
                allele_s_p = allele_s2
                blockid_s_p = blockid_s2 
                changes_varpos_next.append(varpos)

            else:
                print("blockIds of s1 and s2 don't match: ",varpos,allele_s_p,blockid_s_p,)
            
    if allele_s1 != '.' :
        if allele_s1 != allele_s2 and allele_s2 == allele_s3 and allele_s3 == allele_s4:
            result_consensus.append(allele_s1 == allele_s5)
            #print(varpos, allele_s1 == allele_s5)
            allele_s_p = str(1-int(allele_s_p))
        
    if allele_s1 == '.' and allele_s5 != ".":
    #if allele_s5 != ".":
        
        
        cntr0+=1
        if allele_s2 != allele_s3 and allele_s3 == allele_s4 and allele_s4 == allele_s5:
            cntr1+=1
        if allele_s2 == allele_s3 and allele_s3 == allele_s4 and allele_s4 == allele_s5:
            cntr2+=1
        if allele_s2 == allele_s3 and allele_s3 == allele_s4 and allele_s4 != allele_s5:
            cntr3+=1
        if allele_s2 != allele_s3 and allele_s3 == allele_s4 and allele_s4 != allele_s5:
            cntr4+=1
        
        
print("s5", cntr0)         
print("s2!=s3=s4=s5:",cntr1,"s2=s3=s4=s5:",cntr2)
print("s2=s3=s4!=s5:",cntr3,"s2!=s3=s4!s5:",cntr4) 
        
#     matrix_row_list_prd = [str(varpos),str(allele_s1)+":"+str(blockid_s1),
#                                   str(allele_s2)+":"+str(blockid_s2),str(allele_s3)+":"+str(blockid_s3),
#                                   str(allele_s4)+":"+str(blockid_s4), str(allele_s5),
#                                   str(allele_s_p)+":"+str(blockid_s_p)] 
    
#     matrix_row_lists_prd.append(matrix_row_list_prd)
    
    
    


# In[ ]:


len(matrix_row_lists)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


varpos_list_inmatrix=[]
for matrix_row_list in matrix_row_lists_prd:
    varpos_list_inmatrix.append(matrix_row_list[0])
print(len(varposs_all),len(varpos_list_inmatrix))


# In[ ]:





# In[ ]:


var_lines_list_imp = []
for var_line_list in var_lines_list: # 
    var_line_list_part=var_line_list[:9]
    
        
#     matrix_row_list_prd = [str(varpos),str(allele_s1)+":"+str(blockid_s1),
#                                   str(allele_s2)+":"+str(blockid_s2),str(allele_s3)+":"+str(blockid_s3),
#                                   str(allele_s4)+":"+str(blockid_s4), str(allele_s5),
#                                   str(allele_s_p)+":"+str(blockid_s_p)] 
    varpos = var_line_list_part[1]
    if varpos in varpos_list_inmatrix:
        
        gt_flags = var_line_list[8]
        idx_ps = gt_flags.split(":").index("PS")
        idx_gt = gt_flags.split(":").index("GT")        
        s1 = var_line_list[9] # the rest of var info is got from the s1
        split_s1 = s1.split(":")
        
        var_idx = varpos_list_inmatrix.index(varpos)
        matrix_row_list_prd = matrix_row_lists_prd[var_idx]

        allele_s_p = matrix_row_list_prd[6][0]
        blockid_s_p = matrix_row_list_prd[6][2:]
        if allele_s_p != ".":
            split_s1[idx_gt] = allele_s_p+"|"+str(1-int(allele_s_p))
            split_s1[idx_ps] = blockid_s_p 
        else:
            split_s1[idx_gt]="./."
            split_s1[idx_ps]="."
            
        var_line_list_imp = split_s1
        var_line_list_part.append(":".join(var_line_list_imp))
        allele_s5 = matrix_row_list_prd[5]

        if allele_s5 != ".":
        
            var_line_list_part.append(allele_s5+"|"+str(1-int(allele_s5))+":.:.:.:.:.:.:.:.:1")
        else:
            var_line_list_part.append("./.:.:.:.:.:.:.:.:.:.")
            
        var_lines_list_imp.append(var_line_list_part)
        
        


# In[ ]:


var_line_list_imp


# In[ ]:


len(var_line_list_imp)


# In[ ]:


var_lines_list_imp[0]


# In[ ]:


gt_flags


# In[ ]:





# In[ ]:





# In[ ]:


var_line_list[8]


# In[ ]:


matrix_row_list_prd[6][:]


# In[ ]:





# In[ ]:





# In[ ]:


vcf_out_address=vcf_file_address+"_imp.vcf"

vcf_out_file=open(vcf_out_address,'w');

for header_line in header_lines_list:
    vcf_out_file.write(header_line)
    
for var_line_list in var_lines_list_imp:
    
    vcf_out_file.write("\t".join(var_line_list)+"\n")

vcf_out_file.close()
   
    



# In[ ]:


var_line_list_part


# In[ ]:


var_lines_list_hifi = []
for var_line_list in var_lines_list: # 
    var_line_list_part=var_line_list[:9]
    
        
#     matrix_row_list_prd = [str(varpos),str(allele_s1)+":"+str(blockid_s1),
#                                   str(allele_s2)+":"+str(blockid_s2),str(allele_s3)+":"+str(blockid_s3),
#                                   str(allele_s4)+":"+str(blockid_s4), str(allele_s5),
#                                   str(allele_s_p)+":"+str(blockid_s_p)] 
    varpos = var_line_list_part[1]
    if varpos in varpos_list_inmatrix:
        
        gt_flags = var_line_list[8]
        idx_ps = gt_flags.split(":").index("PS")
        idx_gt = gt_flags.split(":").index("GT")        
        s1 = var_line_list[9] # the rest of var info is got from the s1
        split_s1 = s1.split(":")

        blockid_s1 = split_s1[idx_ps]
        allele_s1 = split_s1[idx_gt]
        
        if allele_s1 != "./.":

            var_idx = varpos_list_inmatrix.index(varpos)
            matrix_row_list_prd = matrix_row_lists_prd[var_idx]   

            var_line_list_imp = split_s1
            var_line_list_part.append(":".join(var_line_list_imp))

            allele_s5 = matrix_row_list_prd[5]
            if allele_s5 != ".":
                var_line_list_part.append(allele_s5+"|"+str(1-int(allele_s5))+":.:.:.:.:.:.:.:.:1")
            else:
                var_line_list_part.append("./.:.:.:.:.:.:.:.:.:.")

            var_lines_list_hifi.append(var_line_list_part)
        
        


# In[ ]:


vcf_out_address=vcf_file_address+"_hifi.vcf"

vcf_out_file=open(vcf_out_address,'w');

for header_line in header_lines_list:
    vcf_out_file.write(header_line)
    
for var_line_list in var_lines_list_hifi:
    
    vcf_out_file.write("\t".join(var_line_list)+"\n")

vcf_out_file.close()
   
    



# In[ ]:





# In[ ]:





# In[ ]:



diff_list=[]

for block_i in range(len(allele_blocks_s1)): # 
   
    allele_block_s1=allele_blocks_s1[block_i]
    varpos_block_s1=varpos_blocks_s1[block_i]


    for var_i, varpos in enumerate(varpos_block_s1):
        allele_s1= allele_block_s1[var_i]
        
        if varpos in varposs_s5:
            allele_s5= alleles_s5[varposs_s5.index(varpos)]
        
            if allele_s1!=allele_s5:
                
                for block_i in range(number_phaseblocks_s2): # 
   
                    allele_block_s2=allele_blocks_s2[block_i]
                    varpos_block_s2=varpos_blocks_s2[block_i]
                    if varpos in varpos_block_s2:
                        allele_s2= allele_block_s2[varpos_block_s2.index(varpos)]
                        diff_list.append([varpos,allele_s1,allele_s2,allele_s5])
                        #print([varpos,allele_s1,allele_s2,allele_s5],"\n")
                        break

                
                
                
                
                


# In[ ]:





# In[ ]:



diff_list_s2=[]

for block_i in range(len(allele_blocks_s2)): # 
   
    allele_block_s2=allele_blocks_s2[block_i]
    varpos_block_s2=varpos_blocks_s2[block_i]


    for var_i, varpos in enumerate(varpos_block_s2):
        allele_s2= allele_block_s2[var_i]
        
        if varpos in varposs_s5:
            allele_s5= alleles_s5[varposs_s5.index(varpos)]
        
            if allele_s2!=allele_s5:
                
                for block_i in range(number_phaseblocks_s1): # 
   
                    allele_block_s1=allele_blocks_s1[block_i]
                    varpos_block_s1=varpos_blocks_s1[block_i]
                    if varpos in varpos_block_s1:
                        allele_s1= allele_block_s1[varpos_block_s2.index(varpos)]
                        diff_list_s2.append([varpos,allele_s1,allele_s2,allele_s5])
                        #print([varpos,allele_s1,allele_s2,allele_s5],"\n")
                        break

            


# In[ ]:





# In[ ]:





# In[ ]:


hifiok_illno=[]
hifino_illno=[]
hifino_illok=[]

for [varpos, hifi,illu,parnt] in diff_list_s2:
    if illu!=parnt and hifi==parnt:
        hifiok_illno.append(varpos)
    if illu!=parnt and hifi!=parnt:
        #print([varpos, hifi,illu,parnt] )
        hifino_illno.append(varpos)       
    if illu==parnt and hifi!=parnt:
        #print([varpos, hifi,illu,parnt] )
        hifino_illok.append(varpos) 
    
print(len(diff_list_s2),len(hifiok_illno),len(hifino_illok),len(hifino_illno))


# In[ ]:


hifiok_illno=[]
hifino_illno=[]
hifino_illok=[]

for [varpos, hifi,illu,parnt] in diff_list:
    
    if illu!=parnt and hifi==parnt:
        hifiok_illno.append(varpos)
    if illu!=parnt and hifi!=parnt:
        #print([varpos, hifi,illu,parnt] )
        hifino_illno.append(varpos)       
    if illu==parnt and hifi!=parnt:
        #print([varpos, hifi,illu,parnt] )
        hifino_illok.append(varpos) 
    
print(len(diff_list_s2),len(hifiok_illno),len(hifino_illok),len(hifino_illno))


# In[ ]:


allele_s5==allele_s2
[allele_s1,allele_s5,allele_s2]


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


similiar_count=0
dissimiliar_count=0
dissimiliar_count_s1ok=0
dissimiliar_count_s2ok=0
dissimiliar_count_exist5=0
similiar_count_s1s2ok=0
similiar_count_s1s2no=0
similiar_count_exist5=0

for block_i in range(number_phaseblocks_s2): # 
   
    allele_block_s2=allele_blocks_s2[block_i]
    varpos_block_s2=varpos_blocks_s2[block_i]

    intersected_idins1_s2=intersected_blocks_idins1_s2[block_i]
    if intersected_idins1_s2 !=0:
        block_indx_s1= blockid_blocks_s1.index(intersected_idins1_s2)
    varpos_block_s1=varpos_blocks_s1[block_indx_s1]
    allele_block_s1=allele_blocks_s1[block_indx_s1]
    
    intersect_varpos=set(varpos_block_s1) & set(varpos_block_s2)

    for varpos in intersect_varpos:
        allele_s1= allele_block_s1[varpos_block_s1.index(varpos)]
        allele_s2= allele_block_s2[varpos_block_s2.index(varpos)]
        if allele_s1==allele_s2: 
            similiar_count+=1
            if varpos in varposs_s5:
                allele_s5= alleles_s5[varposs_s5.index(varpos)]
                similiar_count_exist5+=1
                if allele_s1==allele_s5 and allele_s2==allele_s5:
                    similiar_count_s1s2ok+=1
                if allele_s1!=allele_s5 and allele_s2!=allele_s5:
                    similiar_count_s1s2no+=1
        else:
            #print(varpos)
            dissimiliar_count+=1
            if varpos in varposs_s5:
                allele_s5= alleles_s5[varposs_s5.index(varpos)]
                dissimiliar_count_exist5+=1
                if allele_s1==allele_s5 and allele_s2!=allele_s5:
                    dissimiliar_count_s1ok+=1
                if allele_s1!=allele_s5 and allele_s2==allele_s5:
                    dissimiliar_count_s2ok+=1
                    

            
#     if similiar_count>1 and dissimiliar_count>1:
#         print(similiar_count)

print(similiar_count,dissimiliar_count,dissimiliar_count_exist5,dissimiliar_count_s1ok,dissimiliar_count_s2ok)
print(similiar_count_exist5,similiar_count_s1s2ok,similiar_count_s1s2no)


# In[ ]:


len(allele_block_s2)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:



diff_list=[]
count_simi=0
block_i=1
print("blockid ",blockid_blocks_s1[block_i])

allele_block_s1=allele_blocks_s1[block_i]
varpos_block_s1=varpos_blocks_s1[block_i]

print("Pos","HiFi","Illm","parental","\n")
for var_i, varpos in enumerate(varpos_block_s1):
    allele_s1= allele_block_s1[var_i]

    if varpos in varposs_s5:
        allele_s5= alleles_s5[varposs_s5.index(varpos)]

        if allele_s1!=allele_s5:

            for block_i in range(number_phaseblocks_s2): # 

                allele_block_s2=allele_blocks_s2[block_i]
                varpos_block_s2=varpos_blocks_s2[block_i]
                if varpos in varpos_block_s2:
                    allele_s2= allele_block_s2[varpos_block_s2.index(varpos)]
                    diff_list.append([varpos,allele_s1,allele_s2,allele_s5])
                    print(varpos,allele_s1,allele_s2,allele_s5,"\n")
                    break
        if allele_s1==allele_s5:
            count_simi+=1
            


                
print(count_simi)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


sum1=0
for  a in varpos_blocks_s1:
    sum1+=len(a)
print(sum1)


# In[ ]:


sum2=0
for  a in varpos_blocks_s2:
    sum2+=len(a)


# In[ ]:


sum2


# In[ ]:





# In[ ]:





# In[ ]:





# In[2]:


# coverage check

bed_address="/Volumes/work/myjupyter/phaseyou/depth.regions.bed"

bed_file=open(bed_address,'r');

selected_lines=[]
start_pos_list=[]
depth_list=[]

for line in bed_file:
    line_strip=line.strip()
    line_split=line_strip.split("\t")
    
    
    depth=float(line_split[3])
    depth_list.append(depth)
    
    pos=float(line_split[1])
    start_pos_list.append(pos)
    
    
bed_file.close()


# In[8]:


list_pos=[16117006, 16441341, 16442657, 17041656, 18764533, 18770724, 20212954,20234525, 20258973, 20292393, 20321889, 21510471, 25043308, 36007921, 39072681, 48600367, 50428656, 50501957, 51112361, 51146139]

list1_depth=[]
for  pos1 in list_pos:
    depth_pos1=depth_list[int(pos1/500)]
    list1_depth.append(depth_pos1)
    print(pos1/1000000, depth_pos1)
    


# In[9]:


import numpy as np
np.mean(depth_pos1)


# In[10]:


import matplotlib.pyplot as plt  


list_pos1=[i/1e6 for i in list_pos]
plt.bar(list_pos1,list1_depth, color ='maroon', width = .4)
#plt.show()

plt.xlabel("mb")

plt.ylabel("coverage")

plt.savefig("Coverage_wrong_majority.pdf")


# In[ ]:





# In[ ]:




# majority voting




len(matrix_row_lists)


# In[ ]:





# In[ ]:


all_alleles, s1,s2,s3,s4,allele_s5


# In[ ]:


animals = ['cat', 'dog', 'rabbit', 'guinea pig', 'rabbit']

# 'rabbit' is removed
animals.remove('rabbit')

# Updated animals List
print('Updated animals list: ', animals)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


allele_s1='.';blockid_s1='.'
allele_s2='.';blockid_s2='.'
allele_s3='.';blockid_s3='.'
allele_s4='.';blockid_s4='.'
if s1 != '.:.':
    allele_s1 = s1[0]
    blockid_s1 = s1[2:]
    blockid_s1_last_known = blockid_s1
    
if s2 != '.:.':
    allele_s2 = s2[0]
    blockid_s2 = s2[2:]        
if s3 != '.:.':
    allele_s3 = s3[0]
    blockid_s3 = s3[2:]
if s4 != '.:.':
    allele_s4 = s4[0]
    blockid_s4 = s4[2:]

allele_s_p = allele_s1
blockid_s_p = blockid_s1

if allele_s1 == '.' and allele_s2 != ".":

    if blockid_s2 == blockid_s1_last_known:
        allele_s_p = allele_s2
        blockid_s_p = blockid_s2
        changes_varpos_pre.append(varpos)
    else:             # check next variant, 

        blockid_s1_next = matrix_row_lists[row_i+1][1][2:]
        blockid_s2_next = matrix_row_lists[row_i+1][2][2:]
        if blockid_s2_next == blockid_s2:
            allele_s_p = allele_s2
            blockid_s_p = blockid_s2 
            changes_varpos_next.append(varpos)

        else:
            print("blockIds of s1 and s2 don't match: ",varpos,allele_s_p,blockid_s_p,)
        
if allele_s1 != '.' :
    if allele_s1 != allele_s2 and allele_s2 == allele_s3 and allele_s3 == allele_s4:
        result_consensus.append(allele_s1 == allele_s5)
        print(varpos, allele_s1 == allele_s5)
        allele_s_p = str(1-int(allele_s_p))

   
        
        
    
matrix_row_list_prd = [str(varpos),str(allele_s1)+":"+str(blockid_s1),
                              str(allele_s2)+":"+str(blockid_s2),str(allele_s3)+":"+str(blockid_s3),
                              str(allele_s4)+":"+str(blockid_s4), str(allele_s5),
                              str(allele_s_p)+":"+str(blockid_s_p)] 

matrix_row_lists_prd.append(matrix_row_list_prd)
  

