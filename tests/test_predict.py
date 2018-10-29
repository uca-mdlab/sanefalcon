import unittest
import os

class MyTestCase(unittest.TestCase):
    def setUp(self):
        self.scriptname = 'predict.sh'
        self.model = './tests/out.model'
        self.sample = './tests/sample.bam'

        # self.nprev =
        # self.npifwd =
        # self.npirev =
        self.extensions = ['.1','.fwd','.ifwd','.rev','.irev','.cnp','.ff']


    def test1(self):
        nbline2=[]
        outfilecnp = '{}{}'.format(self.sample,self.extensions[5])
        outfileff = '{}{}'.format(self.sample,self.extensions[6])
        result = os.system("bash {} {} {}".format(self.scriptname, self.model,self.sample))
        self.assertEqual(result, 0)

        with open(outfilecnp,'r') as cnpresfile:
            nbline = cnpresfile.readlines()
        cnpresfile.close()
        with open(outfileff,'r') as ffresfile:
            ff = ffresfile.readlines()
        ffresfile.close()
        nbline=[x.split(',\n')[0] for x in nbline]
        ff=[x.split('\n')[0] for x in ff]

        self.assertEqual(len(nbline),2)

        nbline = [x.split(',') for x in nbline]
        if len(nbline) == 2:
            nbelm=nbline[0]+nbline[1]

        self.assertEqual(len(nbelm),294)
        self.assertEqual(len(ff),2)
        self.assertIn('Fetal Fraction',ff[0])


if __name__ == '__main__':
    unittest.main()
