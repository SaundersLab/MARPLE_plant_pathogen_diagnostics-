#!/usr/bin/python
#encoding: utf-8
"""
Created on Tue Nov 20 12:53:43 2012

@author: kyoshida
"""
import sys
import re
class PileUp:
    
    def __init__(self, pileup):
        self.chr, self.position, self.ref, self.depth, self.base, self.quality= pileup.split("\t")
        
    def is_a_snp(self):
        
        base={}
        
        Acount = self.base.count("A")
        Acount += self.base.count("a")
        
        Ccount = self.base.count("C")
        Ccount += self.base.count("c")     
        
        Gcount = self.base.count("G")
        Gcount += self.base.count("g")
        
        Tcount = self.base.count("T")
        Tcount += self.base.count("t")
        
        Refcount = self.base.count(".")
        Refcount += self.base.count(",")
        
        flag=0
        if self.base.find('+') >=0:
            flag = 1
        elif self.base.find('-') >=0:
            flag = 1
        else:
            flag = 0
        
        if flag == 1:
            for m1 in re.finditer('([+|-])(\d+)(\w+)', self.base):
                num = int(m1.group(2))
                b = (m1.group(3))[0:num]          
                b= b.upper()
                Acount -= b.count("A")
                Ccount -= b.count("C")
                Gcount -= b.count("G")
                Tcount -= b.count("T")
        
        base["A"] = Acount
        base["C"] = Ccount
        base["G"] = Gcount
        base["T"] = Tcount
        
        if self.ref == "A":
            base["A"] += Refcount
        elif self.ref == "C":
            base["C"] += Refcount
        elif self.ref == "G":
            base["G"] += Refcount
        elif self.ref == "T":
            base["T"] += Refcount 
        
        
        nucleotide=[]
        depth=[]
        
        for b,d in base.iteritems():
            
            if int(d)>0:
                ratio=float(d)/float(self.depth)
                nucleotide.append(b)
                depth.append(round(ratio,3))
        
        return self.chr, self.position, self.ref, self.depth, ",".join(nucleotide), ",".join(map(str,depth))
        
for pileup in sys.stdin: 
    pileup.rstrip() #remove new line mark

    p = PileUp(pileup)
    (chr,position,ref,depth,nucleotide,ratio) = p.is_a_snp()
        
    print chr+"\t"+position+"\t"+ref+"\t"+depth+"\t"+nucleotide+"\t",ratio;
        
        
        
            
            
            
