import os
import re

def ccs_run(self):
    directory = self.ui.t3le_1.text().strip()
    if os.path.isdir(directory) == False:
        self.ui.t3le_1.setStyleSheet("background-color: rgb(255, 0, 0)")
        return


    #Please do not edit past this point
    files = [x for x in os.listdir(directory) if x.lower().endswith('.mout')]

    opf = open(directory+'//CCS_DATA.csv','w')
    opf.write('Filename,Mass,Temperature,Cross Section,Standard Deviation\n')
    opf.close()

    for file in files:
        opf = open(directory+'//'+file,'r')
        data = opf.read()
        opf.close()
        
        filename = re.findall('input file name = (.*?)\n',data)[0].strip()
        #filename_mod = []
        #for i in filename:
        #    filename_mod+=str(i).strip()

        mass = re.findall('mass of ion = (.*?)\n',data)[0].strip().replace('D','E')

        try:
            CCS_data = ((data[data.find('Mobility Summary')+25:]).replace(' ','').split('\n')[2:])
        except:
            continue

        for data in CCS_data[:-1]:
            data = data.split(',')
            opf = open(directory+'//CCS_DATA.csv','a')
            opf.write('%s,%s,%s,%s,%s\n'%(filename,mass,data[0],data[1],data[2]))
            opf.close()  



        
        
