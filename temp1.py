#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 09:15:27 2018

@author: u1490431
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May 31 16:03:18 2018
@author: James Rig
"""



import numpy as np


sal_d1 = [1,2,3]
sal_d2 = [2,3,4]
sal_d3 = [5,6,7]

pcp_d1 = [2,2,3] 
pcp_d2 = [4,5,8]
pcp_d3 = [5,9,9]


data = [[sal_d1, sal_d2, sal_d3], [pcp_d1, pcp_d2, pcp_d3]]

barscatter(data, transpose=True, paired=True)