# Copyright (C) 2016 VU University Medical Center Amsterdam
# Author: Roy Straver (r.straver@vumc.nl)
#
# This file is part of SANEFALCON
# SANEFALCON is distributed under the following license:
# Attribution-NonCommercial-ShareAlike, CC BY-NC-SA (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
# This license is governed by Dutch law and this license is subject to the exclusive jurisdiction of the courts of the
# Netherlands.


# Arguments:
# 1 Training nucleosome profiles
# 2 Training reference data
# 3 Output basename
# 4 Test nucleosome profiles
# 5 Test reference data

import sys
#import argparse
import glob
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr
import numpy as np
import random
import pickle
import sklearn.linear_model as sklm

# ---------------------------------------------------------------------------- #
# Data loading functions
# ---------------------------------------------------------------------------- #
ignoredSeries = []  # 'serie23_140707','serie12_140428']


def load_nucleosome_file(fname):
    samples = {}
    coverages = {}
    with open(fname) as in_:
        for sample_line in in_:
            arr = sample_line.split(',')
            values = np.array([float(x) for x in arr[1:]])
            sum_values = np.sum(values)
            samples[arr[0]] = np.divide(values, sum_values)
            coverages[arr[0]] = sum_values
    return samples, coverages


def load_reference_file(fname):
    series = dict()
    reference = dict()
    girls = dict()
    bads = dict()
    with open(fname) as in_:
        for line in in_:
            samplename, ffref, gender = line.strip().split("\t")

            if float(ffref) > 30 or float(ffref) < 3:
                continue
            if samplename in ignoredSeries:
                continue
            if gender == 'Male':
                reference[samplename] = float(ffref)
                series[samplename] = "Training"  # splitLine[-1]
            elif gender == 'Female':
                girls[samplename] = float(ffref)
            elif gender == 'BAD':
                bads[samplename] = float(ffref)
            else:
                print(ffref)
    return reference, series, girls, bads


def loadRefFileTrisomy(refFile):
    series = dict()
    reference = dict()
    with open(refFile) as referenceFile:
        for line in referenceFile:
            splitLine = line.strip().split(" ")
            reference[splitLine[1]] = float(splitLine[2])
            series[splitLine[1]] = '21'
    return reference, series


def split_by_reference(samples, reference):
    overlap = [x for x in reference if x in samples]
    overlap.sort()
    noOverlap = [x for x in reference if x not in samples]
    noOverlap.sort()

    print(len(overlap), "samples overlap")
    print(len(noOverlap), "samples noOverlap")
    return overlap, noOverlap


# ---------------------------------------------------------------------------- #
# Data preparation functions
# ---------------------------------------------------------------------------- #
def get_correlation_profile(samples, trainingSet, yVals):
    correlations = []
    for i in range(len(samples[trainingSet[0]])):
        bpVals = [samples[x][i] for x in trainingSet]
        correlations.append(pearsonr(bpVals, yVals)[0])
    return correlations


def get_nucl_ratios(names, sampleSet, correlations):
    selected_regions = []
    sample_scores = []
    for sample_name in names:
        sample = sampleSet[sample_name]
        score = sum([val * correlations[i] for i, val in enumerate(sample)])
        sample_scores.append(score)
        selected_regions.append(sample)
    return sample_scores, selected_regions


def getBinnedProfiles(trainingSet,binSize=25):
    binnedSamples = []
    for sample in trainingSet:
        binValues = []
        for i in range(int(len(sample)/binSize)):
            binValues.append(sum(sample[i * binSize: i * binSize+binSize]))
        binnedSamples.append(binValues)
    return binnedSamples


def get_area_scores(samples):
    scores = []
    for thisSample in samples:
        neg = sum(thisSample[30:60])  # + sum(thisSample[225:250])
        pos = sum(thisSample[80:125])  # + sum(thisSample[160:180])  + sum(thisSample[280:])
        scores.append([pos, neg])
    return scores


# ---------------------------------------------------------------------------- #
# Data processing functions
# ---------------------------------------------------------------------------- #
def get_error_rate(prediction, reference):
    errors = []
    for i, val in enumerate(prediction):
        error = (val - reference[i])
        errors.append(abs(error))
    return np.median(errors)


def train_polyfit(samples, reference):
    return np.polyfit(samples, reference, 1)


def train_linear_model(samples, reference):
    clf = sklm.LinearRegression()  # .Ridge(alpha = .3)#
    clf.fit(samples, reference)
    return clf


def test_polyfit(samples, reference, p, prefix):
    fitSamples=[]
    for i, val in enumerate(samples):
        fitSamples.append(p[0] * val+p[1])
    print(prefix + " polyfit: Pearson:", pearsonr(fitSamples, reference))
    print(prefix + " polyfit: errorRate:", get_error_rate(fitSamples, reference))
    return fitSamples


def test_linear_model(samples, reference, clf, prefix):
    predicted = clf.predict(samples)
    print(prefix + " linearModel: Pearson:", pearsonr(predicted,reference))
    print(prefix + (" linearModel: Residual sum of squares: %.2f" % np.mean((predicted - reference) ** 2)))
    print(prefix + (' linearModel: Variance score: %.2f' % clf.score(samples, reference)))
    print(prefix + " linearModel: errorRate:", get_error_rate(predicted, reference))
    return predicted


def leaveSomeOut(samples, reference, leaveOutSize):
    leaveOutTests = len(samples) / leaveOutSize
    errorRates = []
    for i in range(leaveOutTests):

        tempTrSamp = samples[:i * leaveOutSize]
        tempTrSamp.extend(samples[i * leaveOutSize + leaveOutSize:])
        tempTrRef = reference[:i * leaveOutSize]
        tempTrRef.extend(reference[i * leaveOutSize + leaveOutSize:])

        tempTeSamp = samples[i * leaveOutSize:i * leaveOutSize + leaveOutSize]
        tempTeRef = reference[i * leaveOutSize:i * leaveOutSize + leaveOutSize]

        tempModel = train_linear_model(tempTrSamp, tempTrRef)
        errorRates.append(get_error_rate(tempModel.predict(tempTeSamp), tempTeRef))

    return numpy.mean(errorRates),numpy.std(errorRates)


# ---------------------------------------------------------------------------- #
# Plotting functions
# ---------------------------------------------------------------------------- #
def plotCorrelation(correlations, outFile):
    plt.figure(figsize=(16, 2))
    plt.title("Correlation per BP in Artificial Nucleosome on Training Set")
    plt.ylabel("Pearson-Correlation Score")
    plt.xlabel("Aligned Nucleosome BP Position")

    plt.plot(correlations)
    plt.xlim([0,292])
    center = 147-1
    plt.xticks([0,center-93, center-73, center, center+73, center+93,len(correlations)-1],['\nUpstream','93','73\nStart','0\nCenter','73\nEnd','93','\nDownstream'])

    plt.axvline(x=center-93, linewidth=1, ls='--', color = 'k')
    plt.axvline(x=center-73, linewidth=1, ls='--', color = 'k')
    plt.axvline(x=center, linewidth=1, ls='--', color = 'k')
    plt.axvline(x=center+73, linewidth=1, ls='--', color = 'k')
    plt.axvline(x=center+93, linewidth=1, ls='--', color = 'k')

    plt.savefig(outFile, dpi=100)


def plotScatter(trainX, trainY, testX, testY, outFile, plotName, recolor):

    def fillStdDev(xVals=trainX,yVals=trainY,stdColor='green'):
        # Create a temporary list for sorting data
        tmpList=[]
        for i,val in enumerate(xVals):
            tmpList.append([val,yVals[i]])
        tmpList.sort()

        movAvg=[]
        movStdH=[]
        movStdL=[]
        windowSize=15

        for i, val in enumerate(tmpList):
            local=tmpList[max(0, i-windowSize / 2): i + windowSize / 2 + 1]
            movAvg.append(numpy.mean([x[1] for x in local]))
            movStdH.append(movAvg[-1]+numpy.std([x[1] for x in local]))
            movStdL.append(movAvg[-1]-numpy.std([x[1] for x in local]))
        plt.plot([x[0] for x in tmpList],movAvg,"-",color=stdColor,linewidth="2")
        plt.fill_between([x[0] for x in tmpList[windowSize/2:-windowSize/2]], movStdH[windowSize/2:-windowSize/2], movStdL[windowSize/2:-windowSize/2], color=stdColor, alpha='0.2')

    plt.figure()
    #fillStdDev(trainX,trainY,"green")
    #fillStdDev(trainX,trainY,"blue")
    #fillStdDev(testX,testY,"red")
    recoloredX = [testX[x] for x in recolor]
    recoloredY = [testY[x] for x in recolor]
    normalX = [testX[x] for x in range(len(testX)) if x not in recolor]
    normalY = [testY[x] for x in range(len(testY)) if x not in recolor]
    plt.scatter(trainX,trainY,color="blue")
    #plt.scatter(testX,testY,color="red")
    plt.scatter(normalX,normalY,color="red")
    plt.scatter(recoloredX,recoloredY,color="green")#,marker='s')

    lowerLimit=0#0.0018#0
    upperLimit=25#0.003#25
    #plt.plot([lowerLimit,upperLimit],[lowerLimit,upperLimit],"--",color='black')

    #plt.xlim([lowerLimit,upperLimit])
    #plt.ylim([lowerLimit,upperLimit])
    #plt.xlim([0.00175,0.00325])
    #plt.ylim([0.00175,0.00325])

    plt.xlim([0,25])
    plt.ylim([0,25])


    #plt.xticks([5,10,15,20])
    #plt.yticks([5,10,15,20])

    plt.title(plotName)
    plt.xlabel("FF using Nucleosomes")
    plt.ylabel("FF using Y-Chrom")

    plt.savefig(outFile+".pdf", dpi=100)


def plotProfiles(training,testing,outFile,correlations=[]):
    train=[sum(x) for x in map(list,zip(*training))]
    test=[sum(x) for x in map(list,zip(*testing))]

    train=[x/sum(train) for x in train]
    test=[x/sum(test) for x in test]

    plt.figure(figsize=(16, 2))
    plt.plot(train)
    plt.plot(test,color='red')

    tempCor=[x+1 for x in correlations]
    #plt.plot([x/sum(tempCor) for x in tempCor],color='green')

    plt.xlim([0, 292])
    center = 147-1
    plt.xticks([0,center-93, center-73, center, center+73, center+93,len(train)-1],['\nUpstream','93','73\nStart','0\nCenter','73\nEnd','93','\nDownstream'])

    plt.axvline(x=center-93, linewidth=1, ls='--', color = 'k')
    plt.axvline(x=center-73, linewidth=1, ls='--', color = 'k')
    plt.axvline(x=center, linewidth=1, ls='--', color = 'k')
    plt.axvline(x=center+73, linewidth=1, ls='--', color = 'k')
    plt.axvline(x=center+93, linewidth=1, ls='--', color = 'k')

    plt.title("Aligned Nucleosome Profile")
    plt.xlabel("Aligned Nucleosome BP Position")
    plt.ylabel("Frequency")

    plt.savefig(outFile+".profiles.pdf", dpi=100)


# ---------------------------------------------------------------------------- #
# Main stuff
# ---------------------------------------------------------------------------- #
def main(argv=None):
    if argv is None:
        argv = sys.argv

    print("- Training stage:")
    # Load sample data
    samples, coverages = load_nucleosome_file(argv[1])

    # Load known answers for our data
    reference, series, girls, bads = load_reference_file(argv[2])

    # Filter out samples that are missing in either file
    training_set, training_set_no = split_by_reference(samples, reference)
    covs = [coverages[x] for x in training_set]

    # Match nucleosome profiles with reference values
    yVals = [reference[x] for x in training_set]

    # Obtain the nucleosome correlation profile over the reference samples
    correlations = get_correlation_profile(samples, training_set, yVals)

    smoothedCorrelations = [np.mean(correlations[max(i - 4, 0): i + 5]) for i in range(len(correlations))]
    #correlations=smoothedCorrelations

    # Obtain predictive values
    # sample_scores contains 1 value per sample (correlation weighted sum of read data)
    # selected_regions returns binned values
    sample_scores, profiles = get_nucl_ratios(training_set, samples, correlations)

    selected_regions = get_area_scores(profiles)

    # Fit our models to the data
    poly_fit = train_polyfit(sample_scores, yVals)
    linear_model = train_linear_model(selected_regions, yVals)

    # Test our fit
    tr_polyfit = test_polyfit(sample_scores, yVals, poly_fit, "Train")
    tr_linearmodel = test_linear_model(selected_regions, yVals, linear_model, "Train")

    fittedVals=[0]*len(yVals)
    for i, val in enumerate(sample_scores):
        #print i,val,yVals[i],val*poly_fit[0]+poly_fit[1]
        fittedVals[i] = val * poly_fit[0] + poly_fit[1]
    print(' '.join([str(x) for x in poly_fit]))

    plt.scatter(fittedVals,yVals)
    plt.xlim([0, 25])
    plt.savefig(sys.argv[3]+'.direct2.pdf', dpi=100)

    #quit()

    with open(sys.argv[3]+'.model', 'w') as modelFile:
        modelFile.write(' '.join([str(x) for x in correlations]))
        modelFile.write('\n')
        modelFile.write(' '.join([str(x) for x in poly_fit]))

    # If we were blessed with a test set then use it
    if len(argv) >= 5:
        print ("\n- Testing Stage:")
        # Load test data
        newSamples,newCoverages=loadNuclFile(argv[4])

        # Load known answers for our test data
        newReference,newSeries,newGirls,newBads=loadRefFile(argv[5])

        # Filter out samples that are missing in either file
        testSet,testSetNo = splitByReference(newSamples,newReference)
        newCovs=[newCoverages[x] for x in testSet]

        # Match nucleosome profiles with reference values
        newYVals=[newReference[x] for x in testSet]
        print (testSet)

        singleRun=["15P0008B","15P0135A","15P0136A","15P0137A","15P0138A","15P0139A","15P0140A","15P0148A","15P0149A","15P0150A","15P0151A","15P0152A"]##
        #["15P0184A","15P0185A","15P0186A","15P0187A","15P0197A","15P0198A","15P0199A","15P0200A","15P0201A","15P0206A"]#		print [x for x in testSet if x in singleRun]
        recolor=[]
        for i,val in enumerate(testSet):
            if val in singleRun:
                recolor.append(i)
        #print recolor

        # Obtain predictive values
        # sample_scores contains 1 value per sample (correlation weighted sum of read data)
        # selected_regions returns binned values
        newSampleScores,newSelectedRegions=get_nucl_ratios(testSet, newSamples, correlations)
        newProfiles=newSelectedRegions

        #newSelectedRegions=getBinnedProfiles(newSelectedRegions,bestBinSize)
        newSelectedRegions=get_area_scores(newSelectedRegions)

        # Test our previously trained fits
        tePolyFit     = test_polyfit(newSampleScores, newYVals, poly_fit, "Test")
        teLinearModel = test_linear_model(newSelectedRegions, newYVals, linear_model, "Test")
        plotScatter(tr_polyfit,yVals,tePolyFit,newYVals,argv[3]+".polyfit",'Nucleosome based prediction',recolor)
        plotScatter(tr_linearmodel,yVals,teLinearModel,newYVals,argv[3]+".linearmodel",'linearmodel',recolor)

        #print sample_scores,yVals,newSampleScores,newYVals
        plotScatter(sample_scores,yVals,newSampleScores,newYVals,argv[3]+".direct",'direct',recolor)


        plotScatter(sample_scores,covs,newSampleScores,newCovs,argv[3]+".cov",'coverages',recolor)

        plotProfiles(profiles,newProfiles,argv[3],correlations)

        errorValues=[0]*len(tePolyFit)
        for i in range(len(tePolyFit)):
            #print newYVals[i],tePolyFit[i]
            errorValues[i]=tePolyFit[i]-newYVals[i]
        print (numpy.mean(errorValues),numpy.std(errorValues))
        print (numpy.mean([abs(x) for x in errorValues]),numpy.std([abs(x) for x in errorValues]))


        errorValues=[]
        for i in range(len(tePolyFit)):
            #print newYVals[i],tePolyFit[i]
            if newYVals[i] > 10 and newYVals[i] < 15:
                errorValues.append(tePolyFit[i]-newYVals[i])
        print (numpy.mean(errorValues),numpy.std(errorValues))
        print (numpy.mean([abs(x) for x in errorValues]),numpy.std([abs(x) for x in errorValues]))

        test_polyfit([newSampleScores[x] for x in recolor], [newYVals[x] for x in recolor], poly_fit, "recolor")
        #testLinearModel([newSelectedRegions[x] for x in recolor],[newYVals[x] for x in recolor],linear_model,"recolor")
        plt.figure()
        sampleScoresXX,selectedRegionsXX = get_nucl_ratios(training_set, samples, correlations)
        newSampleScoresXX,newSelectedRegionsXX = get_nucl_ratios(testSetNo, newSamples, correlations)

        #print testPolyFit(sample_scores),testPolyFit(sampleScoresXX),testPolyFit(newSampleScores),testPolyFit(newSampleScoresXX)

        sampleScoresFit	=	[poly_fit[0]*val+poly_fit[1] for val in sample_scores]
        sampleScoresXXFit	=	[poly_fit[0]*val+poly_fit[1] for val in sampleScoresXX]
        newSampleScoresFit	=	[poly_fit[0]*val+poly_fit[1] for val in newSampleScores]
        newSampleScoresXXFit	=	[poly_fit[0]*val+poly_fit[1] for val in newSampleScoresXX]

        #plt.boxplot([sample_scores,sampleScoresXX,newSampleScores,newSampleScoresXX])#,newPredictedNoRef])
        plt.boxplot([sampleScoresFit,sampleScoresXXFit,newSampleScoresFit,newSampleScoresXXFit])#,newPredictedNoRef])
        plt.title("Boxplots for Regressor Fetal Fraction Outputs")
        plt.ylabel("FF score using Nucleosomes")
        plt.xlabel("Dataset used")
        plt.xticks([1,2,3,4], ['Training Boys', 'Training Girls', 'Test Boys','Test Girls'])
        plt.savefig(sys.argv[3]+'.boxplots.pdf', dpi=100)

    # Make some plots to show what we have done
    plotCorrelation(correlations,argv[3]+'.correlations.pdf')
    plotCorrelation(smoothedCorrelations,argv[3]+'.smoothedcorrelations.pdf')
		
if __name__ == "__main__":
	main()

