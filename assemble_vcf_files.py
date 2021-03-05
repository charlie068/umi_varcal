#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 17:37:25 2021

@author: jean-Charles
"""

import os
import glob


rootdir = YourPath

for subdir, dirs, files in os.walk(rootdir):
    print(subdir)
    try:
        
        sample=subdir.split('/')[6]
        print(rootdir+'/'+sample+'_out.vcf')
        with open(rootdir+'/'+sample+'_out.vcf','a+') as fwrite:
            for file in files:
                filename=file.replace('_sorted','')
                ext=filename.split('.')[-1]
                if ext=='vcf':
                    sample=filename.split('.')[0]
                    with open(subdir+'/'+file,'r+') as fread:
                        for line in fread:
                            if line[0]!='#':
                                fwrite.write(line)
                                
    except:
        print('Problem with {}'.format(subdir))
        
