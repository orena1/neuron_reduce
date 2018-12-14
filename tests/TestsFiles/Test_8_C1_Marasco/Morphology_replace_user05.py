name = 'geo5038801mod'
apical_dendriteEnd = 79


total_user5 = 70
f = open(name  + '.hoc','r')

new_ls = ''
for line in f:
    if 'user5[' in line and 'create' not in line and 'append' not in line:
        parts = line.split('user5[')
        #sdfs
        
        
        
        if '{user5[51] connect user5[52](0), 1}' in line:
            pass
            #asdas
            pass
        for i in range(len(parts)):
            if i in range(1,len(parts)):
                #asdas
                parts[i]
                num = int(parts[i][:parts[i].index(']')])
                num = num + apical_dendriteEnd
                new_ls += 'apic[' + str(num) + ']' + parts[i][parts[i].index(']')+1:].replace('apical_dendrite','apic')
            else:
                new_ls += parts[i].replace('apical_dendrite','apic')
    

    elif 'create user5' in line:
        new_ls +='\n'
    
    elif 'user5al.append()' in line:
        new_ls +='        for i=' + str(apical_dendriteEnd) + ', ' + \
                                      str(apical_dendriteEnd + total_user5-1)  +' apic[i] user5al.append()\n'


    elif 'user5[i] all.append()' in line:
        new_ls +='\n'
    elif 'apical_dendrite[i] all.append()' in line:
        new_ls += '        for i=0, '+ str(apical_dendriteEnd + total_user5-1)  +' apic[i] all.append()'
        
    elif 'apical_dendrite' in line and 'create' not in line:
        new_ls += line.replace('apical_dendrite','apic').replace('apical_dendrite','apic').replace('apical_dendrite','apic')
        
    
    elif 'create apical_dendrite' in line:
        parts = line.split('apical_dendrite[')
        new_ls += parts[0] + 'apic[' + str(apical_dendriteEnd+total_user5) +  parts[1][parts[1].index(']'):]
        
    
    elif 'template' in line:
        new_ls += line.replace(name, name + 'Mod')
            
    else:
        new_ls += line








f.close()
f1 = open(name + 'Mod.hoc', 'w')
f1.write(new_ls)
f1.close()
