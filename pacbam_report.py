#! /usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import argparse
from argparse import ArgumentParser

def read_params():
	parser = ArgumentParser(description="")
	inputArg = parser.add_argument('-i', '--input', type=str, required=True,
		help="Specify the input file prefix")
	modeArg = parser.add_argument('-m', '--mode', type=int, required=True,
		help="Specify the mode used")
	outputArg = parser.add_argument('-o', '--output', type=str, required=False,
		help="Specify the output file name (Default input.pdf)")
	strandBiasArg = parser.add_argument('-s', '--strandBias', action='store_true', required=False, 
		help="Plots the strand bias distribution")
	args = parser.parse_args()

	if args.mode < 0 or args.mode > 4:
		raise argparse.ArgumentError(modeArg, "Mode must be between 0 and 4")
        
        if args.strandBias and (args.mode not in [0,1]):
                raise argparse.ArgumentError(strandBiasArg, "Strand Bias is available only in modes 0 and 1")

	return args

def pileupReport(filename, binSize,outFile):
	listSize = 1000000
	with open(filename) as pileup:
		next(pileup) #remove header
		
		cumulativeCoverage = [0] * listSize # << normal experiment coverage
		bins = int(1/binSize) # number of bins for the allelic fraction
		cumulativeAF = [0] * (bins + 1) # allocate bins for the allelic fraction

		totalBases = 0 # bases in pileup
		meanAF = 0 # mean allelic fraction

		for l in pileup: 
			l = l.split('\t') 
			totalBases += 1 # add a base
			coverage = int(l[8]) # 8 is the index of the coverage in the pileup file 
			AF = float(l[7]) # 7 is the index of the Allelic Fraction in the pileup file
			#coverage
			try:
				cumulativeCoverage[coverage] += 1 # add 1 to the element count relative to a specific allelic fraction
			except Exception:
				cumulativeCoverage[len(cumulativeCoverage)-1] += 1 # very high coverage add to the last bin

			#AF

			i = int(AF*(1/binSize)) # index for the allelic fraction (if binsize is multiple of 10 e.g AF = 0.003 and binsize = 0.001 -> i = 3)
			meanAF += AF # add the AF to the mean
			cumulativeAF[i] += 1 # add 1 to the count of the allelic fraction

		#coverage
		cumulativeDistribution = float(cumulativeCoverage[0]) #first element for the cumulative distribution 

		cumulativePointsX = [0] # array to be plotted (x axis) coverage value
		cumulativePointsY = [1] # array to be plotted (y axis) cumulative distribution
		xl = 0
		for i in range(1,listSize): # size of the list 
			cumulativeDistribution += cumulativeCoverage[i] # add the element to the cumulative distribution
			cumulativePointsX.append(i) #add the coverage value 
			cumulativePointsY.append(1 - cumulativeDistribution/totalBases) #add the cumulative distribution up to point i
			if 1 - (cumulativeDistribution/totalBases) < 0.01 and xl == 0:
				xl = i

		#make the plot
		plt.plot(cumulativePointsX,cumulativePointsY)
		plt.title("Cumulative Coverage Distribution")
		plt.xlim(0,xl)
		plt.ylim(0,1)
		plt.ylabel("Base at coverage")
		plt.xlabel("Coverage")
		plt.text(.0, -.1, "File: " + filename,transform=plt.gca().transAxes)
		plt.savefig(outFile,format='pdf')
		plt.clf()
		
		#Plot AF
		zeroAF = cumulativeAF[0]
		meanAF = meanAF/totalBases
		plt.hist([x * binSize for x in range(1,bins+1)],bins=bins,weights=cumulativeAF[1:len(cumulativeAF)])
		plt.title("Allelic Fraction Distribution")
		plt.text(.0, -.1, "File: " + filename,transform=plt.gca().transAxes)
		plt.text(0.5,0.9,s=("Total Bases: %d\n Bases with reference Allele: %d\n Bases with reference Allele (fraction): %.4f\n Mean Allelic Fraction: %.4f") % (totalBases,zeroAF,float(zeroAF)/totalBases,meanAF), horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes, bbox=dict(alpha=0))
		plt.xlabel("AF")
		plt.ylabel("Read count")
		plt.xlim(0,1)
		
		plt.savefig(outFile,format='pdf')
		plt.clf()

		#Plot AF zoom up to 0.2
		plt.hist([x * binSize for x in range(1,bins+1)],bins=bins,weights=cumulativeAF[1:len(cumulativeAF)])
		plt.title("Allelic Fraction Distribution (Zoom: 0.2)")
		plt.text(.0, -.1, "File: " + filename,transform=plt.gca().transAxes)
		plt.text(0.5,0.9,s=("Total Bases: %d\n Bases with reference Allele: %d\n Bases with reference Allele (fraction): %.4f\n Mean Allelic Fraction: %.4f") % (totalBases,zeroAF,float(zeroAF)/totalBases,meanAF), horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes, bbox=dict(alpha=0))
		plt.xlabel("AF")
		plt.ylabel("Read count")
		plt.xlim(0,0.2)
		plt.savefig(outFile,format='pdf')
		plt.clf()
		
def SNPsReport(filename,lowerBound,upperBound,outFile):

	with open(filename) as snps:
		totalBins = [0,0,0] #3 Bins
		next(snps) #skip header 

		coverageDistribution = []

		for l in snps:
			l = l.split('\t') 
			AF = float(l[9]) #get allelic fraction
			coverage = int(l[10])
			coverageDistribution.append(coverage)
			#divide the value in the 3 totalBins 0 <= AF < lowerbound, lowerbound < AF <= upperBound, AF > upperbound
			if AF <= lowerBound:
				totalBins[0] += 1
			elif AF <= upperBound:
				totalBins[1] += 1
			else:
				totalBins[2] += 1

		# 3 Bins (for labels)
		q1 = (np.percentile(coverageDistribution,25))
		q2 = (np.percentile(coverageDistribution,50))
		q3 = (np.percentile(coverageDistribution,75))
		q4 = (np.percentile(coverageDistribution,100))

		with open(filename) as snps:
			stratBins = [[0,0],[0,0],[0,0],[0,0]] #bins for each base for heterozygous and homozygous
			next(snps) #skip header 

			for l in snps:
				l = l.split('\t') 
				AF = float(l[9]) #get allelic fraction
				coverage = int(l[10])
				#check the coverage first quintile second third or forth
				if coverage < q1:
					binsIndex = 0
				elif coverage <= q2:
					binsIndex = 1
				elif coverage <= q3:
					binsIndex = 2
				elif coverage <= q4:
					binsIndex = 3

				# check the allelic fraction for heterozygous or homozygous
				if AF <= upperBound and AF > lowerBound:
					stratBins[binsIndex][0] += 1
				elif AF > upperBound:
					stratBins[binsIndex][1] += 1

		#barplot of SNPs
		pET = plt.bar(0,totalBins[1], color=EColor) 
		pOT = plt.bar(0,totalBins[2],bottom=totalBins[1], color=OColor) 

		pEQ1 = plt.bar(1,stratBins[0][0], color=EColor) 
		POQ1 = plt.bar(1,stratBins[0][1],bottom=stratBins[0][0], color=OColor)

		pEQ2 = plt.bar(2,stratBins[1][0], color=EColor) 
		POQ2 = plt.bar(2,stratBins[1][1],bottom=stratBins[1][0], color=OColor)

		pEQ3 = plt.bar(3,stratBins[2][0], color=EColor) 
		POQ3 = plt.bar(3,stratBins[2][1],bottom=stratBins[2][0], color=OColor)

		pEQ4 = plt.bar(4,stratBins[3][0], color=EColor) 
		POQ4 = plt.bar(4,stratBins[3][1],bottom=stratBins[3][0], color=OColor)

		Epatch = mpatches.Patch(color=EColor, label='Alternative Heterozygous')
		Opatch = mpatches.Patch(color=OColor, label='Alternative Homozygous')
		plt.text(.0, -.134, "File: " + filename + "\nReference Homozygous: AF $\leq$ " + str(lowerBound) + "\nAlternative Heterozygous: " + str(lowerBound) + " < AF $\leq$ " + str(upperBound) +  "\nAlternative Homozygous: AF > " + str(upperBound) ,transform=plt.gca().transAxes)
		plt.text(0.5,0.95,s=("Total SNPs: %d\nReference Homozygous SNPs: %d") % (sum(totalBins),totalBins[0]),horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes, bbox=dict(alpha=0))
		plt.xticks(range(5), ["Total", "CoverageQ1", "CoverageQ2", "CoverageQ3", "CoverageQ4"])
		plt.ylabel("Count")

		plt.legend(handles=[Epatch,Opatch],bbox_to_anchor=(1, 1))
		plt.title("SNPs types")
		plt.savefig(outFile,format='pdf')
		plt.clf()


def SNVsReport(filename,outFile,strandBias):

	with open(filename) as snvs:
		next(snvs) #skip header 
		#list of bases (for index)
		bases = ["A","C","G","T"]
		#possible changes A = 0 C = 1 G = 2 T = 3
		positionCombination = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]] 

		strandBiasDistribution = [] #initialize the strand bias distribution
		
		for l in snvs:
			l = l.split('\t') 
			
			refBase = l[2] # get the reference base
			altBase = l[3] # get the alternative base
			
			#ACGT order
			baseCount = [0,0,0,0] #get the base count 
			baseCount[0] = int(l[4])
			baseCount[1] = int(l[5])
			baseCount[2] = int(l[6])
			baseCount[3] = int(l[7])		

			refIndex = bases.index(refBase) # get the index of the reference base

			for b in bases: #for each base update the counter for the change (if present)
				if b != refBase:
					if b == "A" and baseCount[0] != 0:
						positionCombination[refIndex][0] += 1 #baseCount[0]
					elif b == "C" and baseCount[1] != 0:
						positionCombination[refIndex][1] += 1 #baseCount[1]
					elif b == "G" and baseCount[2] != 0:
						positionCombination[refIndex][2] += 1 #baseCount[2]
					elif b == "T" and baseCount[3] != 0:
						positionCombination[refIndex][3] += 1 #baseCount[3]

			if strandBias:
			#//////if strandBias option////////
				#compute the strandBias
				reverseBaseCount = [0,0,0,0]
				reverseBaseCount[0] = int(l[10])
				reverseBaseCount[1] = int(l[11])
				reverseBaseCount[2] = int(l[12])
				reverseBaseCount[3] = int(l[13])


				if altBase != 'N' and baseCount[bases.index(refBase)] > 0:
					refReverse = float(reverseBaseCount[bases.index(refBase)])
					refForward = float(baseCount[bases.index(refBase)] - refReverse)
					 
					altReverse = float(reverseBaseCount[bases.index(altBase)])
					altForward = float(baseCount[bases.index(altBase)] - altReverse)

					#strandBias Score
					strandBiasValue = abs((refReverse/baseCount[bases.index(refBase)])-(altReverse/baseCount[bases.index(altBase)]))
					#add the score to the distribution
					strandBiasDistribution.append(strandBiasValue)
				
		if strandBias:
			#strandbias plot
			plt.hist(strandBiasDistribution,bins=[x/100.0 for x in range(0,100)])
			plt.title("Strand Bias Distribution")
			plt.axvline(np.percentile(strandBiasDistribution,25), color=Q1Color, linestyle='dashed', linewidth=2)
			plt.axvline(np.percentile(strandBiasDistribution,50), color=medianColor, linestyle='dashed', linewidth=2)
			plt.axvline(np.percentile(strandBiasDistribution,75), color=Q3Color, linestyle='dashed', linewidth=2)
			Q1Patch = mlines.Line2D([],[],color=Q1Color, label='Q1',linestyle='dashed')
			medianPatch = mlines.Line2D([],[],color=medianColor, label='Median',linestyle='dashed')
			Q3Patch = mlines.Line2D([],[],color=Q3Color, label='Q3',linestyle='dashed')
			plt.legend(handles=[Q1Patch,medianPatch,Q3Patch],bbox_to_anchor=(1, 1))
			plt.xlabel("Strand bias")
			plt.ylabel("Count")
			plt.text(.0, -.1, "File: " + filename,transform=plt.gca().transAxes)
			plt.savefig(outFile,format='pdf')
			plt.clf()


		#base modification plot
		pAC = plt.bar(0,positionCombination[0][1], color=CColor) # c
		pAG = plt.bar(0,positionCombination[0][2],bottom=positionCombination[0][1], color=GColor) # g
		pAT= plt.bar(0,positionCombination[0][3],bottom=positionCombination[0][1]+positionCombination[0][2],color=TColor) 		

		pCA = plt.bar(1,positionCombination[1][0],color=AColor)
		pCG = plt.bar(1,positionCombination[1][2],bottom=positionCombination[1][0],color=GColor)
		pCT= plt.bar(1,positionCombination[1][3],bottom=positionCombination[1][0]+positionCombination[1][2],color=TColor)

		pGA = plt.bar(2,positionCombination[2][0],color=AColor)
		pGC = plt.bar(2,positionCombination[2][1],bottom=positionCombination[2][0],color=CColor)
		pGT= plt.bar(2,positionCombination[2][3],bottom=positionCombination[2][0]+positionCombination[2][1],color=TColor)

		pTA = plt.bar(3,positionCombination[3][0],color=AColor)
		pTC = plt.bar(3,positionCombination[3][1],bottom=positionCombination[3][0],color=CColor)
		pTG= plt.bar(3,positionCombination[3][2],bottom=positionCombination[3][0]+positionCombination[3][1],color=GColor)

		Apatch = mpatches.Patch(color=AColor, label='A')
		Cpatch = mpatches.Patch(color=CColor, label='C')
		Gpatch = mpatches.Patch(color=GColor, label='G')
		Tpatch = mpatches.Patch(color=TColor, label='T')
		plt.xticks(range(4), bases)
		plt.title("Bases modification types")
		plt.ylabel("Count")
		plt.xlabel("Reference Base")
		plt.legend(handles=[Apatch,Cpatch,Gpatch,Tpatch],bbox_to_anchor=(1, 1),title="Alternative")
		plt.text(.0, -.1, "File: " + filename,transform=plt.gca().transAxes)

		plt.savefig(outFile,format='pdf')
		plt.clf()

def RCReport(filename,outFile):

	#////REGION WITH SIZE >5///////

	with open(filename) as rc:
		next(rc) #skip header

		regionCoverageGlobal = []
		rCoverage = []
		gcDistribution = []
		for l in rc:
			l = l.split('\t') #split
			#considering regions with size at least 6
			if float(l[2]) - float(l[1]) > 5: #////REGION WITH SIZE >5///////
				gcDistribution.append(float(l[7]))
				regionCoverageGlobal.append(float(l[5]))
				rCoverage.append(float(l[6]))

		#gc content plot
		plt.hist(gcDistribution,bins=[x/100.0 for x in range(0,100)])
		plt.axvline(np.percentile(gcDistribution,25), color=Q1Color, linestyle='dashed', linewidth=2)
		plt.axvline(np.percentile(gcDistribution,50), color=medianColor, linestyle='dashed', linewidth=2)
		plt.axvline(np.percentile(gcDistribution,75), color=Q3Color, linestyle='dashed', linewidth=2)		
		Q1Patch = mlines.Line2D([],[],color=Q1Color, label='Q1',linestyle='dashed')
		medianPatch = mlines.Line2D([],[],color=medianColor, label='Median',linestyle='dashed')
		Q3Patch = mlines.Line2D([],[],color=Q3Color, label='Q3',linestyle='dashed')
		plt.legend(handles=[Q1Patch,medianPatch,Q3Patch],bbox_to_anchor=(1, 1))
		plt.title("GC content")
		plt.xlabel("GC content")
		plt.ylabel("Count")
		plt.text(.0, -.1, "File: " + filename,transform=plt.gca().transAxes)
		plt.savefig(outFile,format='pdf')
		plt.clf()

		
		#region coverage plot		
		f, (ax1, ax2) = plt.subplots(1,2,sharey=True,sharex=True)

		ax1.set_ylabel("Region at coverage")
		ax1.set_xlabel("Coverage")
		ax2.set_xlabel("Coverage")
		
		ax1.hist(regionCoverageGlobal, bins=[x for x in range(0,int(np.percentile(regionCoverageGlobal,99)))],cumulative=-1,density=True)
		ax1.set_title("Region Coverage Global")
		ax2.hist(rCoverage, bins= [x for x in range(0,int(np.percentile(regionCoverageGlobal,99)))], cumulative=-1,density=True)
		ax2.set_title("Region Coverage Peaked")

		ax1.text(-1.25, -.1, "File: " + filename,transform=plt.gca().transAxes)
		
		plt.suptitle("Region Coverage Distribution")
		plt.savefig(outFile,format='pdf')
		plt.clf()




#color constant
OColor = "0.25"
EColor = "0.75"

Q1Color = "0.75"
medianColor = "0.5"
Q3Color = "0.25"

AColor = "#006DDB"
CColor = "#490092"
GColor = "#920000"
TColor = "#B6DBFF"


#////////////////MAIN//////////////////////
if __name__ == '__main__':

	plt.rcParams['figure.figsize'] = 10.5, 8.7 # size of figures
	args = read_params()

	if args.output is not None:
		pp = PdfPages(args.output) # outfile
	else:
		pp = PdfPages(args.input + ".pdf")

	#input files
	pileupFile = args.input + ".pileup" 
	snpsFile = args.input + ".snps"
	snvsFile = args.input + ".snvs"
	RCFile = args.input + ".rc"

	#mode
	mode = args.mode

	#report on the pilup
	if mode == 1 or mode == 4:
		print "Computing Pileup Report"
		pileupReport(pileupFile, 0.001,pp)

	#report on snp
	if mode == 0 or mode == 1 or mode == 2:
		print("Computing SNPs Report")
		SNPsReport(snpsFile,0.2,0.8,pp)

	#report on SNVs
	if mode == 0 or mode == 1:
		print "Computing SNVs Report"
		SNVsReport(snvsFile,pp,args.strandBias)

	#report rc
	if mode == 0 or mode == 1 or mode == 3:
		print "Computing RC Report"
		RCReport(RCFile,pp)

	pp.close()#close the output file
