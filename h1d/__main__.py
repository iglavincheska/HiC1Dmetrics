import argparse
import pandas as pd
from .MultiTypeScore import *
from .plotMetrics import *
from .plotTwoSample import *
from .MultiSampleScore import *
from .callDirectionalTAD import *
from .calldTADAllchr import *
from .discrete import *
from .stripe import *
import os

def CLI():
    parser = argparse.ArgumentParser(description="HiC1Dmetrics is python3-based tools to \
                                                calculate, visualize, and analyze 1D metrics for Hi-C data \n \
                                                (https://github.com/wangjk321/HiC1Dmetrics) ")
    #parser.set_defaults(func=lambda args: parser.print_help())
    subparsers = parser.add_subparsers(help="Choose the mode to use sub-commands")

    #Function 1
    def func_basic(args):
        if not os.path.exists(args.data):print("File does not exist");exit(1)
        args.matrix = args.data
        if args.mode == "plot":
            if args.datatype == "rawhic":
                path = hic2matrix(args.matrix,args.resolution,args.chromosome,args.gt,args.juicertool)
                if args.controlmatrix:
                    controlpath = hic2matrix(args.controlmatrix,args.resolution,args.chromosome,args.gt,args.juicertool)
                else:
                    controlpath = None
            else:
                path = args.matrix
                controlpath = args.controlmatrix

            if args.plottype == "tri":
                if not controlpath:
                    PlotTri(path,args.resolution,args.chromosome,startSite=args.start,endSite=args.end).draw()
                elif controlpath:
                    DiffDraw(path,controlpath,args.resolution,args.chromosome,startSite=args.start,endSite=args.end).draw_tri()
            elif args.plottype == "square":
                if not controlpath:
                    PlotSquare(path,args.resolution,args.chromosome,startSite=args.start,endSite=args.end).draw()
                elif controlpath:
                    DiffDraw(path,controlpath,args.resolution,args.chromosome,startSite=args.start,endSite=args.end).draw_square()
            elif args.plottype == "tad":
                if not controlpath:
                    PlotTAD(path,args.resolution,args.chromosome,startSite=args.start,endSite=args.end).drawTAD(squareSize=int(args.parameter))
                elif controlpath:
                    DiffDraw(path,controlpath,args.resolution,args.chromosome,startSite=args.start,endSite=args.end).drawTAD(squareSize=int(args.parameter))
            plt.savefig(args.outname+".pdf")
        elif args.mode == "dump":

            codepath = os.path.dirname(os.path.realpath(__file__))
            makeIntra = codepath+"/extract/makeMatrixIntra.sh"
            if not args.juicertool:
                juicer = codepath+"/jc/jctool_1.11.04.jar"
            else:
                juicer = args.juicertool

            if args.chromosome == "all":
                allJuicer(args.data,args.normalize,args.resolution,args.gt,args.outname,juicer,
                        maxchr=args.maxchr,num=args.nProcesser)
                exit(0)

            if args.datatype not in ["rawhic"] or not args.gt:
                print("Error: dump requires rawhic file and genome_table file"); exit(1)



            foldername = args.outname
            os.system("bash "+makeIntra+" "+args.normalize+" "+"."+" "+args.matrix+" "+
                    str(args.resolution)+" "+args.gt+" "+juicer+" "+args.chromosome+" "+foldername)
        elif args.mode == "gd":
            codepath = os.path.dirname(os.path.realpath(__file__))
            gdcode = codepath+"/gd/makeDensity.sh"
            os.system("bash "+gdcode+" -r "+str(args.resolution)+" -g "+args.data+" -t "+args.chromosome+" -o "+args.outname)
        else:
            print("Unsupported mode"); exit(1)
    parser_basic = subparsers.add_parser("basic",help="Provide basic functions to visualize and handle Hi-C data.")
    parser_basic.add_argument('mode', type=str, help='Type of 1D metrics,,should be one of {dump,plot}')
    parser_basic.add_argument('data', type=str, help='Path of matrix file from JuicerResult')
    parser_basic.add_argument('resolution', type=int,help="Resolution of input matrix")
    parser_basic.add_argument("chromosome",type=str,help="Chromosome number.")
    parser_basic.add_argument("-o","--outname",help="output name",type=str,default="unname")
    parser_basic.add_argument('-c','--controlmatrix', type=str, help='Path of control matrix file from JuicerResult',default=None)
    parser_basic.add_argument('--datatype',type=str,help="matrix or rawhic",default="matrix")
    parser_basic.add_argument('--gt',type=str,help="genome table",default="")
    parser_basic.add_argument('--plottype',type=str,help="Type of plot, could be one of {tri,square,tad}",default="tri")
    parser_basic.add_argument('-s','--start',type=int,help="Start sites for plotting",default=0)
    parser_basic.add_argument('-e','--end',type=int,help="End sites for plotting",default=0)
    parser_basic.add_argument("-p","--parameter",type=str,help="Parameter for indicated metrics",default=None)
    parser_basic.add_argument("--normalize",type=str,help="Normalize methods {NONE/VC/VC_SQRT/KR}",default="KR")
    parser_basic.add_argument("-n","--nProcesser",type=int,help="Number of processors",default=10)
    parser_basic.add_argument('--maxchr',type=int,help="Maximum index of chromosome (human genome is 22,i.e.)",default=None)
    parser_basic.add_argument('--juicertool',type=str,help="Specify juicertool with different version.",default=None)
    parser_basic.set_defaults(func=func_basic)

    #Function 2
    #=============================================================================
    def func_one(args):
        if not os.path.exists(args.data):print("File does not exist");exit(1)
        
        # IG: Auto-convert .cool/.mcool to .hic for IF calculation
        if args.type == "IF" and (args.data.endswith('.cool') or args.data.endswith('.mcool')):
            from .loadfile import cool2hic
            print("Detected .cool file for IF calculation, converting to .hic...")
            if not args.gt:
                print("Error: genome_table (-gt) is required for .cool to .hic conversion"); exit(1)
            args.data = cool2hic(args.data, args.resolution, args.gt, args.juicertool)
            args.datatype = "rawhic"
            print(f"Conversion complete. Using .hic file: {args.data}")
        
        args.matrix = args.data
        if args.parameter and args.type not in ["PC1","IF"]:
            args.parameter = int(args.parameter)

        if args.chromosome == "all":
            if not args.maxchr: print("Please sepcify the maximum chromosome"); exit(1)
            if not os.path.exists(args.data):
                print("path not exist"); exit(1)
            if args.datatype == "matrix":
                scoreAll = oneScoreAllchr(args.data,args.resolution,args.type,args.parameter,
                                          maxchr=args.maxchr,prefix=args.prefix,num=args.nProcesser,juicer=args.juicertool)
            else:
                # For rawhic/cool input, compute each chromosome from the same source file.
                if not args.gt:
                    print("Error: --gt is required when chromosome=all and datatype is not matrix"); exit(1)
                gt_df = pd.read_csv(args.gt, sep="\t", header=None, dtype={0: str})
                chrom_list = [str(c) for c in gt_df.iloc[:, 0].tolist()
                              if str(c).lower() not in ["mt", "m", "chrm", "chrmt"]]
                chrom_list = chrom_list[:args.maxchr]
                if len(chrom_list) == 0:
                    print("Error: no chromosomes found in genome table"); exit(1)

                all_scores = []
                for chrom in chrom_list:
                    score = multiScore(args.matrix,args.resolution,chrom,juicer=args.juicertool
                                       ).obtainOneScore(args.type,parameter=args.parameter,datatype=args.datatype,
                                                        gt=args.gt,TADfile=args.TADfile,msi=args.msi)
                    all_scores.append(score)
                scoreAll = pd.concat(all_scores, axis=0, ignore_index=True)
            scoreAll.to_csv(args.outname+"_"+args.type+"_allchr.csv",sep="\t",index=False)
            scoreAll.to_csv(args.outname+"_"+args.type+"_allchr.bedGraph",sep="\t",header=False,index=False)
            exit(0)

        score = multiScore(args.matrix,args.resolution,args.chromosome,juicer=args.juicertool
                            ).obtainOneScore(args.type,parameter=args.parameter,datatype=args.datatype,
                            gt=args.gt,TADfile=args.TADfile,msi=args.msi)
        print("Saving...")
        score.to_csv(args.outname + ".bedGraph", sep="\t", header=False, index=False)

        if args.draw:
            if args.chromosome=="all": print("Error: not supported"); exit(1)
            if args.type == "IF":
                args.parameter = args.outname + ".bedGraph"
            print("==========output figure==========")
            PlotBedGraph(args.matrix,args.resolution,args.chromosome,startSite=args.start,
                        endSite=args.end,datatype=args.datatype,
                        gt=args.gt,juicer=args.juicertool).draw(args.type,UniqueParameter=args.parameter)
            plt.savefig(args.outname+".pdf")


    parser_one = subparsers.add_parser("one",help="1D metrics designed for one Hi-C sample.",
                                            description="1D metrics designed for one Hi-C sample.")
    parser_one.add_argument('type', type=str, help='Type of 1D metrics,,should be one of {IS,CI,DI,SS,DLR,PC1,IES,IAS,IF}.')
    parser_one.add_argument('data', type=str, help='Path of matrix or rawhic file.')
    parser_one.add_argument('resolution', type=int,help="Resolution of input matrix.")
    parser_one.add_argument("chromosome",type=str,help="Chromosome number.")
    parser_one.add_argument("-p","--parameter",type=str,help="Parameter for indicated metrics.",default=None)
    parser_one.add_argument("-o","--outname",help="output name (default: 'metrics').",type=str,default="unname")
    parser_one.add_argument("-d","--draw",action='store_true',help="Plot figure for candidate region.",default=False)
    parser_one.add_argument("-s",'--start',type=int,help="Start sites for plotting.",default=0)
    parser_one.add_argument("-e",'--end',type=int,help="End sites for plotting.",default=0)
    parser_one.add_argument('--datatype',type=str,help="Type of input data: [matrix(default),rawhic,cool].",default="matrix")
    parser_one.add_argument('--msi',type=str,help="Method for significant interactions: [fithic2,hiccups]",default="fithic2")
    parser_one.add_argument('--gt',type=str,help="genome_table file.",default="")
    parser_one.add_argument('--prefix',type=str,help="${prefix}chr1.matrix.gz",default="observed.KR.")
    parser_one.add_argument('--maxchr',type=int,help="Maximum index of chromosome (human genome is 22,i.e.)",default=None)
    parser_one.add_argument("-n","--nProcesser",type=int,help="Number of processors",default=10)
    parser_one.add_argument("-t","--TADfile",type=str,help="Give a TAD file, instead of using building-in TAD calling method",default=None)
    parser_one.add_argument('--juicertool',type=str,help="Specify juicertool with different version.",default=None)
    parser_one.set_defaults(func=func_one)

    #Function 3
    #=============================================================================
    def func_two(args):
        if not os.path.exists(args.data):print("File does not exist");exit(1)
        args.matrix = args.data
        args.controlmatrix = args.controldata

        if args.chromosome == "all":
            if not args.maxchr: print("Please sepcify the maximum chromosome"); exit(1)
            if not os.path.exists(args.data):
                print("path not exist"); exit(1)
            scoreAll =twoScoreAllchr(args.data,args.controldata,args.resolution,args.type,args.parameter,
                                    maxchr=args.maxchr,prefix=args.prefix)
            scoreAll.to_csv(args.outname+"_"+args.type+"_allchr.csv",sep="\t",index=False)
            exit(0)

        ms = multiScore(args.matrix,args.resolution,args.chromosome,control_path=args.controlmatrix,juicer=args.juicertool)
        score,path,control_path = ms.obtainTwoScore(args.type,parameter=args.parameter,datatype=args.datatype,gt=args.gt)
        print("Saving...")
        score.to_csv(args.outname + ".bedGraph", sep="\t", header=False, index=False)

        if args.draw:
            if args.type != "IFC":
                DiffDraw(path,control_path,args.resolution,chr = args.chromosome,startSite=args.start,
                        endSite=args.end,datatype="matrix",gt=args.gt).drawMetric("custom",customfile=args.outname + ".bedGraph",name=args.type)
            elif args.type == "IFC":
                DiffDraw(path,control_path,args.resolution,chr = args.chromosome,startSite=args.start,
                        endSite=args.end,datatype="rawhic",gt=args.gt).drawMetric("custom",customfile=args.outname + ".bedGraph",name=args.type)
            plt.savefig(args.outname+".pdf")

        os.system("rm -rf MatrixTemp*")

    parser_two = subparsers.add_parser("two",help="1D metrics designed for comparison of two Hi-C samples",
                                            description="1D metrics designed for comparison of two Hi-C samples")
    parser_two.add_argument('type', type=str, help='Type of 1D metrics for two-sample comparison,should be one of {ISC,CIC,SSC,deltaDLR,CD,IESC,IASC,IFC,DRF}')
    parser_two.add_argument('data', type=str, help='Path of treated file (matrix or rawhic).')
    parser_two.add_argument('controldata', type=str, help='Path of control file (matrix or rawhic).')
    parser_two.add_argument('resolution', type=int,help="Resolution of input matrix")
    parser_two.add_argument("chromosome",type=str,help="Chromosome number ('chr21',i.e).")
    parser_two.add_argument("-p","--parameter",type=str,help="Parameter for indicated metrics",default=None)
    parser_two.add_argument("-o","--outname",help="output name (default: metricsChange)",type=str,default="metricsChange")
    parser_two.add_argument("-d","--draw",action='store_true',help="Plot figure for candidate region",default=False)
    parser_two.add_argument('-s','--start',type=int,help="Start sites for plotting",default=0)
    parser_two.add_argument('-e','--end',type=int,help="End sites for plotting",default=0)
    parser_two.add_argument('--datatype',type=str,help="Type of input data: [matrix(default),rawhic,cool].",default="matrix")
    parser_two.add_argument('--gt',type=str,help="genome table file",default="")
    parser_two.add_argument('--prefix',type=str,help="${prefix}chr1.matrix.gz",default="observed.KR.")
    parser_two.add_argument('--maxchr',type=int,help="Maximum index of chromosome (human genome is 22,i.e.)",default=None)
    parser_two.add_argument("-n","--nProcesser",type=int,help="Number of processors",default=10)
    parser_two.add_argument('--juicertool',type=str,help="Specify juicertool with different version.",default=None)
    parser_two.set_defaults(func=func_two)

    #Function 4
    #=============================================================================
    def func_types(args):
        if not os.path.exists(args.data):print("File does not exist");exit(1)
        args.matrix = args.data
        typelist = args.typelist.split(",")
        parameterlist = args.parameter.split(",")
        if not args.controlmatrix:
            if not set(typelist).issubset(["IS","CI","DI","SS","DLR","PC1","IES","IAS","IF"]):
                print("Error: not supported"); exit(1)
            if "IF" in typelist and args.datatype == "matrix":
                print("Error: IF required rawhic or cool datatype"); exit(1)
            if "IF" in typelist and args.datatype == "cool":
                if not args.gt:
                    print("Error: genome_table is required for cool to hic conversion"); exit(1)
                print("Detected IF in multitypes with cool input, converting to .hic...")
                args.matrix = cool2hic(args.matrix, args.resolution, args.gt, args.juicertool)
                args.datatype = "rawhic"
            ms = multiScore(args.matrix,args.resolution,args.chromosome,juicer=args.juicertool)
            if not args.draw:
                score = ms.allOneScore(typelist,parameterlist,datatype=args.datatype,gt=args.gt)
            elif args.draw:
                score = ms.plotOneScore(typelist,parameterlist,datatype=args.datatype,gt=args.gt,start=args.start,end=args.end)
                plt.savefig(args.outname+".pdf")
            score.to_csv(args.outname + ".csv", sep="\t", header=True, index=False)
        elif args.controlmatrix:
            if not set(typelist).issubset(["ISC","CIC","SSC","deltaDLR","CD","IESC","IASC","IFC","DRF"]):
                print("Error: not supported"); exit(1)
            if "IFC" in typelist and args.datatype == "matrix": print("Error: IFC required rawhic or cool datatype"); exit(1)
            if "IFC" in typelist and args.datatype == "cool":
                if not args.gt:
                    print("Error: genome_table is required for cool to hic conversion"); exit(1)
                print("Detected IFC in multitypes with cool input, converting treated/control to .hic...")
                args.matrix = cool2hic(args.matrix, args.resolution, args.gt, args.juicertool)
                args.controlmatrix = cool2hic(args.controlmatrix, args.resolution, args.gt, args.juicertool)
                args.datatype = "rawhic"
            if "DRF" in typelist:
                DRFpos = typelist.index("DRF")
                parameterlist[DRFpos] = parameterlist[DRFpos].split("-")
                print(parameterlist)
            ms = multiScore(args.matrix,args.resolution,args.chromosome,control_path=args.controlmatrix,juicer=args.juicertool)
            if not args.draw:
                score = ms.allTwoScore(typelist,parameterlist,datatype=args.datatype,gt=args.gt)
            elif args.draw:
                score = ms.plotTwoScore(typelist,parameterlist,datatype=args.datatype,gt=args.gt,start=args.start,end=args.end)
                plt.savefig(args.outname+".pdf")
            score.to_csv(args.outname + ".csv", sep="\t", header=True, index=False)

        os.system("rm -rf MatrixTemp*")

    parser_types = subparsers.add_parser("multitypes",help="Various types of 1D metrics for the same sample",
                                            description="Various types of 1D metrics for the same sample")
    parser_types.add_argument('typelist', type=str, help='Type of 1D metrics,should be {IS,CI,DI,SS,DLR,PC1,IES,IAS,IF} or {ISC,CIC,SSC,deltaDLR,CD,IESC,IASC,IFC,DRF}')
    parser_types.add_argument('data', type=str, help='Path of matrix file or raw .hic file')
    parser_types.add_argument('resolution', type=int,help="Resolution of input matrix")
    parser_types.add_argument("chromosome",type=str,help="Chromosome number.")
    parser_types.add_argument("-p","--parameter",type=str,help="Parameter for indicated metrics",default=None,required=True)
    parser_types.add_argument('-c','--controlmatrix', type=str, help='Path of control matrix file from JuicerResult',default=None)
    parser_types.add_argument("-o","--outname",help="output name (default metrics)",type=str,default="multitypes_metrics")
    parser_types.add_argument('--datatype',type=str,help="Type of input data: [matrix(default),rawhic,cool].",default="matrix")
    parser_types.add_argument('--gt',type=str,help="genome table",default="")
    parser_types.add_argument("-d","--draw",action='store_true',help="Plot figure for candidate region",default=False)
    parser_types.add_argument('-s','--start',type=int,help="Start sites for plotting",default=0)
    parser_types.add_argument('-e','--end',type=int,help="End sites for plotting",default=0)
    parser_types.add_argument('--juicertool',type=str,help="Specify juicertool with different version.",default=None)
    parser_types.set_defaults(func=func_types)

    #Function 5
    #=============================================================================
    def func_samples(args):
        if not os.path.exists(args.data):print("File does not exist");exit(1)
        if args.type == "IF" and args.datatype == "matrix": print("Error: IF required rawhic datatype"); exit(1)
        datafile = pd.read_csv(args.data,sep="\t",header=None)
        labels = list(datafile.iloc[:,0])
        samplelist = list(datafile.iloc[:,1])
        if not args.corr and not args.heat and not args.line and not args.discrete and not args.anova:
            score = getMultiSamplesScore(samplelist,labels,args.resolution,args.chromosome,args.type,args.parameter,
                                        datatype=args.datatype,gt=args.gt,TADfile=args.TADfile)
        elif args.corr:
            ms = repQC(samplelist,labels,args.resolution,args.chromosome,args.type,args.parameter,
                        datatype=args.datatype,gt=args.gt,TADfile=args.TADfile)
            score = ms.score
            ms.corr_plot()
            plt.savefig(args.outname+"_corr.pdf")
        elif args.discrete:
            if args.datatype == "rawhic":
                print("not supported");exit(1)
            msd = multiSampleDiscrete(samplelist,labels,args.resolution,args.chromosome,args.type,args.parameter)
            plotpath= samplelist[0]
            msd.plotMultiDiscrete(plotpath,args.start,args.end)
            plt.savefig(args.outname+"_discrete.pdf")
            exit(0)
        elif args.anova:
            if args.end == 0:
                print("Please specify a region to calculate p-value");exit(1)
            ms = repQC(samplelist,labels,args.resolution,args.chromosome,args.type,args.parameter,datatype=args.datatype,gt=args.gt)
            score = ms.score
            plist,qlist = ms.anova_like(args.start,args.end)
            df = pd.DataFrame()
            df["chr"] = [args.chromosome]* len(range(args.start,args.end,args.resolution))
            df["start"] = list(range(args.start,args.end,args.resolution))
            df["end"] = df["start"] + args.resolution
            df['pvalue'] = plist
            df['qvalue'] = qlist
            df.to_csv(args.outname+"_anova.txt",sep="\t",index=False)

        elif args.line or args.heat:
            ms = repQC(samplelist,labels,args.resolution,args.chromosome,args.type,args.parameter,
                        datatype=args.datatype,gt=args.gt,TADfile=args.TADfile)
            score = ms.score
            if args.datatype == "matrix":
                plotpath= samplelist[0]
            elif args.datatype == "rawhic":
                plotpath = hic2matrix(samplelist[0],args.resolution,args.chromosome,args.gt,args.juicertool)

            if args.heat: plottype = "heat"
            elif args.line: plottype = "line"

            ms.heatmap_tri(plotpath,args.start,args.end,clmax=args.clmax,heatmin=None,plottype=plottype)
            plt.savefig(args.outname+"_"+plottype+".pdf")

        score.to_csv(args.outname + ".csv", sep="\t", header=True, index=False)

        #print(score.iloc[550:650,:])
        os.system("rm -rf MatrixTemp*")

    parser_samples = subparsers.add_parser("multisamples",help="The same metrics for muliple samples",
                                            description="The same metrics for muliple samples")
    parser_samples.add_argument('type', type=str, help='Type of 1D metrics,,should be one of {IS,CI,DI,SS,DLR,PC1,IES,IAS,IF} or {ISC,CIC,SSC,deltaDLR,CD,IESC,IASC,IFC,DRF}')
    parser_samples.add_argument('data',type=str,help="a txt contain paths for all samples")
    parser_samples.add_argument('resolution', type=int,help="Resolution of input matrix")
    parser_samples.add_argument("chromosome",type=str,help="Chromosome number.")
    parser_samples.add_argument('--datatype',type=str,help="matrix or rawhic",default="matrix")
    parser_samples.add_argument('--samplelist', type=str, help='list of file path, can be rawhic or matrix')
    parser_samples.add_argument('--labels', type=str, help='list of file name')
    parser_samples.add_argument("-p","--parameter",type=str,help="Parameter for indicated metrics",default=None)
    parser_samples.add_argument("-o","--outname",help="output name (default metrics)",type=str,default="multisamples_metrics")
    parser_samples.add_argument('--gt',type=str,help="genome table",default="")
    parser_samples.add_argument("--corr",action='store_true',help="Plot correlation for all samples",default=False)
    parser_samples.add_argument("--heat",action='store_true',help="Plot raw heatmap for all samples",default=False)
    parser_samples.add_argument("--line",action='store_true',help="Plot line chart for all samples",default=False)
    parser_samples.add_argument("--anova",action='store_true',help="Plot line chart for all samples",default=False)
    parser_samples.add_argument("--discrete",action='store_true',help="Plot discrete heatmap for all samples",default=False)
    parser_samples.add_argument('-s','--start',type=int,help="Start sites for plotting",default=0)
    parser_samples.add_argument('-e','--end',type=int,help="End sites for plotting",default=0)
    parser_samples.add_argument('--clmax',type=int,help="End sites for plotting",default=None)
    parser_samples.add_argument("-t","--TADfile",type=str,help="Give a TAD file, instead of using building-in TAD calling method",default=None)
    parser_samples.add_argument('--juicertool',type=str,help="Specify juicertool with different version.",default=None)

    parser_samples.set_defaults(func=func_samples)

    #Function 6
    #=============================================================================
    def func_call(args):
        if not os.path.exists(args.data):print("File does not exist");exit(1)
        args.matrix = args.data
        if args.datatype in ["rawhic",'cool'] and args.mode not in ["hubs", "IFregions"]:
            if args.datatype == "rawhic":
                path = hic2matrix(args.matrix,args.resolution,args.chromosome,args.gt,args.juicertool)
                if args.controlmatrix: controlpath = hic2matrix(args.controlmatrix,args.resolution,args.chromosome,args.gt)
            elif args.datatype == "cool":
                path = cool2matrix(args.matrix,args.resolution,args.chromosome,args.gt)
                if args.controlmatrix: controlpath = cool2matrix(args.controlmatrix,args.resolution,args.chromosome,args.gt)
        else:
            path = args.matrix
            controlpath = args.controlmatrix

        if args.mode == "TAD":
            if not args.parameter: args.parameter = 300000
            tad = TADcallIS(path,args.resolution,args.chromosome,squareSize=int(args.parameter))
            tad.to_csv(args.outname + "_TAD.csv", sep="\t", header=True, index=False)
        elif args.mode == "dTAD":
            if args.parameter:
                parameter = args.parameter.split("-")
            else:
                parameter = [200000,5000000]
            if not args.controlmatrix: print("Error: DRF requires control data"); exit(1)
            dt = DirectionalTAD(path,controlpath,args.resolution,args.chromosome,
                                startDRF=int(parameter[0]), sizeDRF=int(parameter[1]))
            leftdTAD,rightdTAD,_ = dt.extractRegion()
            leftdTAD.to_csv(args.outname + "_leftdTAD.csv", sep="\t", header=True, index=False)
            rightdTAD.to_csv(args.outname + "_rightdTAD.csv", sep="\t", header=True, index=False)
        elif args.mode == "stripeTAD":
            st = stripeTAD(path,args.resolution,args.chromosome)
            sTAD = st.callStripe(squareSize=300000)
            sTAD.to_csv(args.outname + "_stripeTAD.csv", sep="\t", header=True, index=False)
        elif args.mode == "stripe":
            st = call_stripe(path,args.resolution,args.chromosome)
            stripe = st.callStripe(squareSize=300000,strong_thresh=0.2)
            stripe.to_csv(args.outname + "_stripes.csv", sep="\t", header=True, index=False)
        elif args.mode == "hubs":
            if args.datatype != "rawhic": print("Error: hubs requires rawhic datatype"); exit(1)
            if args.chromosome == "all":
                if not args.gt:
                    print("Error: hubs with chromosome=all requires genome table (--gt)"); exit(1)
                gt_df = pd.read_csv(args.gt, sep="\t", header=None, dtype={0: str})
                chrom_list = [str(c) for c in gt_df.iloc[:, 0].tolist()
                              if str(c).lower() not in ["mt", "m", "chrm", "chrmt"]]
                if len(chrom_list) == 0:
                    print("Error: no chromosomes found in genome table"); exit(1)

                if_list = []
                for chrom in chrom_list:
                    score = InteractionFrequency(path,args.resolution,chrom,gt=args.gt,juicer=args.juicertool).getIF()
                    if_list.append(score)
                IF = pd.concat(if_list, axis=0, ignore_index=True)

                # Call hubs per chromosome to avoid large-chromosome bias.
                thresh = IF.groupby("chr").InteractionFreq.transform(lambda s: np.nanpercentile(s, 90))
                hubregion = IF[IF.InteractionFreq > thresh]
            else:
                IF = InteractionFrequency(path,args.resolution,args.chromosome,gt=args.gt,juicer=args.juicertool).getIF()
                thresh = np.percentile(IF.iloc[:,3],90)
                hubregion = IF[IF.iloc[:,3]>thresh]

            hubregion.to_csv(args.outname + "_hubs_IF.csv", sep="\t", header=True, index=False)
            os.system("sed '1d' " + args.outname + "_hubs_IF.csv" + " | sort -k1,1 -k2,2n | bedtools merge -i stdin > "+ args.outname + "_hubs.csv")
        elif args.mode == "IFregions":
            # Read normalized IF bedGraph and extract high-IF regions.
            if not os.path.exists(args.data):
                print("File does not exist"); exit(1)

            threshold = float(args.parameter) if args.parameter else 1.5
            merge_gap = max(0, int(args.resolution) - 1)

            score = pd.read_csv(args.data, sep="\t", header=None, comment="#", dtype={0: str})
            if score.shape[1] < 4:
                print("Error: IFregions expects at least 4 columns: chr,start,end,IFscore"); exit(1)

            score = score.iloc[:, :4].copy()
            score.columns = ["chr", "start", "end", "IF"]
            score["chr"] = score["chr"].astype(str)
            score["start"] = pd.to_numeric(score["start"], errors="coerce")
            score["end"] = pd.to_numeric(score["end"], errors="coerce")
            score["IF"] = pd.to_numeric(score["IF"], errors="coerce")
            score = score.dropna(subset=["start", "end", "IF"])

            high = score[score["IF"] > threshold].copy()
            high.to_csv(args.outname + "_IFregions_raw.bedGraph", sep="\t", header=False, index=False)

            if high.shape[0] == 0:
                pd.DataFrame(columns=["chr", "start", "end"]).to_csv(
                    args.outname + "_IFregions_merged.bed", sep="\t", header=False, index=False
                )
            else:
                high = high.sort_values(["chr", "start", "end"])
                merged = []
                cur_chr = None
                cur_start = None
                cur_end = None

                for _, row in high.iterrows():
                    r_chr = str(row["chr"])
                    r_start = int(row["start"])
                    r_end = int(row["end"])

                    if cur_chr is None:
                        cur_chr, cur_start, cur_end = r_chr, r_start, r_end
                        continue

                    if r_chr == cur_chr and (r_start - cur_end) <= merge_gap:
                        cur_end = max(cur_end, r_end)
                    else:
                        merged.append([cur_chr, cur_start, cur_end])
                        cur_chr, cur_start, cur_end = r_chr, r_start, r_end

                merged.append([cur_chr, cur_start, cur_end])
                pd.DataFrame(merged, columns=["chr", "start", "end"]).to_csv(
                    args.outname + "_IFregions_merged.bed", sep="\t", header=False, index=False
                )
        else:
            print("unsupported model");exit(1)

    parser_call = subparsers.add_parser("call",help="Extract secondary information from metrics (dTAD, stripeTAD, et.al)",
                                            description="Extract secondary information from metrics (dTAD, stripeTAD, et.al)")
    parser_call.add_argument('mode', type=str, help='Running mode,,should be one of {dTAD,stripe,stripeTAD,TAD,hubs,IFregions}')
    parser_call.add_argument('data', type=str, help='Path of matrix file or raw .hic file')
    parser_call.add_argument('resolution', type=int,help="Resolution of input matrix")
    parser_call.add_argument("chromosome",type=str,help="Chromosome number.")
    parser_call.add_argument("-o","--outname",help="output name",type=str,default="defaultname")
    parser_call.add_argument('-c','--controlmatrix', type=str, help='Path of control matrix file from JuicerResult',default=None)
    parser_call.add_argument('--datatype',type=str,help="Type of input data: [matrix(default),rawhic,cool].",default="matrix")
    parser_call.add_argument('--gt',type=str,help="genome table",default="")
    parser_call.add_argument("-p","--parameter",type=str,help="Parameter for indicated mode",default=None)
    parser_call.add_argument('--juicertool',type=str,help="Specify juicertool with different version.",default=None)
    parser_call.set_defaults(func=func_call)

    #Function 7
    #=============================================================================
    def _split_csv_field(value):
        if value is None:
            return []
        text = str(value).strip()
        if text == "" or text.lower() == "nan":
            return []
        return [x.strip() for x in text.split(",") if x.strip()]

    def _to_int_or_none(value):
        if value is None:
            return None
        text = str(value).strip()
        if text == "" or text.lower() == "nan":
            return None
        return int(float(text))

    def _resolve_juicer_path(juicer_tool_value, juicer_dir_value):
        tool = (juicer_tool_value or "").strip() if juicer_tool_value else ""
        if tool:
            return tool

        jdir = (juicer_dir_value or "").strip() if juicer_dir_value else ""
        if not jdir:
            return None

        if os.path.isfile(jdir):
            return jdir

        if not os.path.isdir(jdir):
            print("Error: juicerdir is not a directory or jar file -> " + jdir)
            exit(1)

        preferred = [
            os.path.join(jdir, "juicer_tools.2.20.00.jar"),
            os.path.join(jdir, "jctool_1.11.04.jar"),
        ]
        for p in preferred:
            if os.path.exists(p):
                return p

        jars = [f for f in os.listdir(jdir) if f.endswith(".jar")]
        if len(jars) > 0:
            return os.path.join(jdir, jars[0])

        print("Error: no Juicer jar found in juicerdir -> " + jdir)
        exit(1)

    def func_batch(args):
        if not os.path.exists(args.config):
            print("File does not exist"); exit(1)

        table = pd.read_csv(args.config, sep="\t", dtype=str).fillna("")
        if table.shape[0] == 0:
            print("Error: empty batch file"); exit(1)

        required = ["data", "resolution"]
        for col in required:
            if col not in table.columns:
                print("Error: batch file is missing required column: " + col)
                exit(1)

        for idx, row in table.iterrows():
            data = row.get("data", "").strip()
            if data == "" or not os.path.exists(data):
                print(f"Error (row {idx+1}): invalid data path -> {data}")
                exit(1)

            res = _to_int_or_none(row.get("resolution", ""))
            if res is None:
                print(f"Error (row {idx+1}): invalid resolution")
                exit(1)

            dtype = row.get("datatype", "").strip() or args.datatype
            ext = os.path.splitext(data)[1].lower()
            # Be permissive for batch inputs: if path is .hic, treat it as rawhic even if datatype says cool.
            if ext == ".hic" and dtype == "cool":
                print(f"[batch] row {idx+1}: datatype=cool with .hic input detected, switching datatype to rawhic")
                dtype = "rawhic"
            if ext in [".cool", ".mcool"] and dtype == "rawhic":
                print(f"[batch] row {idx+1}: datatype=rawhic with cool input detected, switching datatype to cool")
                dtype = "cool"
            gt = row.get("gt", "").strip() or args.gt
            juicer = _resolve_juicer_path(
                row.get("juicertool", "").strip() or args.juicertool,
                row.get("juicerdir", "").strip() or row.get("juicer_dir", "").strip() or args.juicerdir,
            )
            tadfile = row.get("TADfile", "").strip() or row.get("tadfile", "").strip()
            genedensity = row.get("genedensity", "").strip() or row.get("geneDensity", "").strip()
            msi = row.get("msi", "").strip() or args.msi

            type_field = row.get("type", "").strip() or row.get("types", "").strip() or row.get("track", "").strip()
            if type_field == "":
                print(f"Error (row {idx+1}): missing type/types column value")
                exit(1)
            typelist = _split_csv_field(type_field)

            param_field = row.get("parameter", "").strip() or row.get("parameters", "").strip()
            params = _split_csv_field(param_field)

            chrom_field = row.get("chromosome", "").strip()
            maxchr = _to_int_or_none(row.get("maxchr", ""))
            if maxchr is None:
                maxchr = args.maxchr

            outprefix = row.get("outname", "").strip() or row.get("outprefix", "").strip()
            if outprefix == "":
                outprefix = f"batch_row{idx+1}"

            # Build chromosome list.
            # If chromosome is empty or 'all', use chromosome names from genome table to preserve naming convention.
            if chrom_field == "" or chrom_field.lower() == "all":
                if not gt:
                    print(f"Error (row {idx+1}): chromosome is empty/all but no genome table (gt) was provided")
                    exit(1)
                gt_df = pd.read_csv(gt, sep="\t", header=None, dtype={0: str})
                chrom_list = [str(c) for c in gt_df.iloc[:, 0].tolist()
                              if str(c).lower() not in ["mt", "m", "chrm", "chrmt"]]
                if maxchr:
                    chrom_list = chrom_list[:maxchr]
            else:
                chrom_list = _split_csv_field(chrom_field)
                if len(chrom_list) == 0:
                    chrom_list = [chrom_field]

            for t_i, mode in enumerate(typelist):
                mode_param = params[t_i] if t_i < len(params) else ""
                if mode == "PC1" and mode_param == "" and genedensity:
                    mode_param = genedensity

                print(f"[batch] row {idx+1}: start metric {mode}")

                calc_dtype = dtype
                calc_data = data

                # Auto-convert cool to hic for IF in batch mode.
                if mode == "IF" and calc_dtype == "cool":
                    if not gt:
                        print(f"Error (row {idx+1}): IF with cool input requires gt")
                        exit(1)
                    in_ext = os.path.splitext(calc_data)[1].lower()
                    if in_ext == ".hic":
                        calc_dtype = "rawhic"
                    elif in_ext in [".cool", ".mcool"]:
                        calc_data = cool2hic(calc_data, res, gt, juicer)
                        calc_dtype = "rawhic"
                    else:
                        print(f"Error (row {idx+1}): IF with datatype=cool requires .cool/.mcool/.hic input, got {calc_data}")
                        exit(1)

                all_scores = []
                for chrom in chrom_list:
                    parameter = mode_param if mode_param != "" else None
                    if parameter and mode not in ["PC1", "IF"]:
                        parameter = int(parameter)

                    score = multiScore(calc_data, res, chrom, juicer=juicer).obtainOneScore(
                        mode,
                        parameter=parameter,
                        datatype=calc_dtype,
                        gt=gt,
                        TADfile=tadfile if tadfile else None,
                        msi=msi,
                    )
                    all_scores.append(score)

                if len(all_scores) == 1:
                    out_df = all_scores[0]
                    out_path = f"{outprefix}_{mode}.bedGraph"
                    out_df.to_csv(out_path, sep="\t", header=False, index=False)
                    print(f"[batch] row {idx+1}: finished metric {mode} -> {out_path}")
                else:
                    out_df = pd.concat(all_scores, axis=0, ignore_index=True)
                    out_csv = f"{outprefix}_{mode}_allchr.csv"
                    out_bg = f"{outprefix}_{mode}_allchr.bedGraph"
                    out_df.to_csv(out_csv, sep="\t", index=False)
                    out_df.to_csv(out_bg, sep="\t", header=False, index=False)
                    print(f"[batch] row {idx+1}: finished metric {mode} -> {out_bg}")

            print(f"Finished batch row {idx+1}")

    parser_batch = subparsers.add_parser("batch", help="Run batch one-sample jobs from a TSV config file")
    parser_batch.add_argument("config", type=str, help="Tab-separated batch config file")
    parser_batch.add_argument("--datatype", type=str, help="Default datatype when missing in rows", default="matrix")
    parser_batch.add_argument("--gt", type=str, help="Default genome table when missing in rows", default="")
    parser_batch.add_argument("--juicertool", type=str, help="Default juicer tools jar when missing in rows", default=None)
    parser_batch.add_argument("--juicerdir", type=str, help="Default Juicer directory (or jar path) when missing in rows", default=None)
    parser_batch.add_argument("--msi", type=str, help="Default IF method [fithic2,hiccups]", default="fithic2")
    parser_batch.add_argument("--maxchr", type=int, help="Optional limit for chromosomes taken from genome table", default=None)
    parser_batch.set_defaults(func=func_batch)

    parser.add_argument("-V","--version",help="Show h1d version",action='store_true',default=False)
    args = parser.parse_args()
    if args.version:
        print("h1d version 0.2.9")
        exit(0)
    try:
        func = args.func
    except AttributeError:
        parser.error("too few arguments")
    func(args)
    os.system("rm -rf MatrixTemp*")
    os.system("rm -rf info.txt")

if __name__ == '__main__':
    CLI()
