# iSANEFALCON

### Improved Single reAds Nucleosome-basEd FetAL fraCtiON

This software is forked from [Sanefalcon](https://github.com/VUmcCGP/sanefalcon). Please refer to the [original README](SANEFALCON.md) for further information.

#### Motivation 
  * Automatization
  * Performance 
  * Parallelization
  * Python3.6+
  
#### Usage

Edit `sanefalcon.conf` with appropriate paths and locations where you `bam` files are stored.

*Optional*: edit lines `64-65` in `manager/file_manager.py` to order the samples in a desired manner.
The goal is to group together samples sequenced on the same microarray (run) in the same batch, to avoid biases.  

Then, prepare a `filename.txt` in the form:

```
samplename1
samplename2
samplename3
...
```
`filaneme.txt` contains the names of the samples for which you want to calculate the fetal fraction.

Then, run:

`python3 main.py filename.txt`