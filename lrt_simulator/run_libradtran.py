from string import Template 
import subprocess as subproc
import io
import os

class uvspec:

    def __init__(self):

        self.exe="../../libRadtran-2.0.4/bin/uvspec"
        self.template_file='uvspec_template.txt'
        self.data_path="../../libRadtran-2.0.4/examples"
        self.io_dir="lrt_io/"
            
        self.options={}    
        self.options["sza"]=0.0
        self.options["umu"]=1.0
        self.options["phi"]=0.0 
        self.options["zout"]=1.0
        self.options["fluorescence"]="fluorescence_file "+self.data_path+"/UVSPEC_FLUORESCENCE.FLU"
        #self.options["fluorescence"]="fluorescence 0.0"
        self.options["albedo"]="albedo_file "+self.data_path+"/UVSPEC_FLUORESCENCE.TOC"
        #self.options["albedo"]="albedo 1.0"

    def gen_input_filename(self):
        """placeholder
        """
        fn='./uvspec'
        for opt in self.options:
            tmp=str(self.options[opt]).replace("../","")
            tmp=tmp.replace("/","_DIR_")
            tmp=tmp.replace(".","p")
            tmp=tmp.replace(" ","_")
            fn+="--"+opt+"_"+tmp   
        return self.io_dir+fn+'.inp'

    def write_uvspec_input(self,overwrite=False):
        """write a new uvspecinput file
        based on the template file and the
        options in the class"""

        inpfilename=self.gen_input_filename()
        if os.path.isfile(inpfilename) and not overwrite:
            return inpfilename

        with open(self.template_file, 'r') as f:
            uvspec_template = f.read()
        with open(inpfilename, 'w') as f:
            print(Template(uvspec_template).safe_substitute(self.options),file=f)

        return inpfilename

    def run(self,overwrite=False):
        
        inpfilename=self.write_uvspec_input(overwrite=overwrite)
        #trick to replace the last occurrence of "inp" with "out":
        outfilename=inpfilename[::-1].replace("pni","tuo",1)[::-1]
        if os.path.isfile(outfilename) and not overwrite:
            return outfilename

        finp=open(inpfilename, 'r')
        fout=open(outfilename, 'w')
        p=subproc.call([self.exe],stdin=finp,stdout=fout)

        return outfilename


if __name__=="__main__":

    u=uvspec()
    u.options["sza"]=3.14777
    #print(u.gen_input_filename())
    #u.write_uvspec_input()
    u.run(overwrite=False)
    
    
