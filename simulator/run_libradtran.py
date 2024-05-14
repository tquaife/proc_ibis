from string import Template 
import subprocess as subproc
import io


class uvspec:

    def __init__(self):

        self.exe="../../libRadtran-2.0.4/bin/uvspec"
        self.template_file='uvspec_template.inp'
        self.data_path="../../libRadtran-2.0.4/examples"
            
        self.options={}    
        self.options["sza"]=30.
        self.options["umu"]=1.0
        self.options["phi"]=0.0 
        self.options["zout"]=1.0
        #self.options["fluorescence"]="fluorescence_file "+data_path+"/UVSPEC_FLUORESCENCE.FLU"
        self.options["fluorescence"]="fluorescence 0.0"
        #self.options["albedo"]="albedo_file "+data_path+"/UVSPEC_FLUORESCENCE.TOC"
        self.options["albedo"]="albedo 1.0"

    def gen_input_filename(self):
        """placeholder
        """
        return'uvspec_TEST.inp'

    def write_uvspec_input(self):
        """write a new uvspecinput file
        based on the template file and the
        options in the class"""

        inpfilename=self.gen_input_filename()
        with open(self.template_file, 'r') as f:
            uvspec_template = f.read()
        with open(inpfilename, 'w') as f:
            print(Template(uvspec_template).safe_substitute(self.options),file=f)

        return inpfilename

    def run(self):
        
        inpfilename=self.write_uvspec_input()
        #trick to replace the last occurrence in "inp" with "out":
        outfilename=inpfilename[::-1].replace("pni","tuo",1)[::-1]
        finp=open(inpfilename, 'r')
        fout=open(outfilename, 'w')
        p=subproc.call([self.exe],stdin=finp,stdout=fout)


if __name__=="__main__":

    u=uvspec()
    u.options["sza"]=3.14777

    u.run()
