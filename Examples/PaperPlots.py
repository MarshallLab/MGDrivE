import MoNeT_MGDrivE as monet

#############################################################################
colors = ['#9f00cc', '#ec0b43', '#232ed1', '#ff009d', '#000000']
style = {
    "width": .01, "alpha": .05, "dpi": 2*1024, "legend": True,
    "aspect": 20, "colors": colors
}
styleS = {
    "width": 0, "alpha": .85, "dpi": 2*1024, "legend": True,
    "aspect": .175, "colors": colors
}
cmaps = monet.generateAlphaColorMapFromColorArray(colors)
#############################################################################
basePath = "/Users/sanchez.hmsc/Desktop/SoftwarePaper/"
experimentString = "ReplacementStochastic"
aggregationDictionary = monet.generateAggregationDictionary(
    ["W", "H", "R", "B"],
    [
        [0, 0, 1, 2, 3],
        [1, 4, 4, 5, 6],
        [2, 5, 7, 7, 8],
        [3, 6, 8, 9, 9]
    ]
)
#############################################################################
repsFolders = monet.listDirectoriesWithPathWithinAPath(
    basePath + experimentString + "/"
)
landscapeReps = monet.loadAndAggregateLandscapeDataRepetitions(
    repsFolders, aggregationDictionary,
    male=False, female=True
)
tracesPlot = monet.plotAllTraces(landscapeReps, style)
monet.quickSaveFigure(
    tracesPlot,
    basePath + experimentString + "_Trace.png", format="png"
)
#############################################################################
filenames = monet.readExperimentFilenames(
    basePath + experimentString + "_Analyzed/"
)
landscapeSumData = monet.sumLandscapePopulationsFromFiles(
    filenames,
    male=True, female=True
)
aggData = monet.aggregateGenotypesInNode(
    landscapeSumData,
    aggregationDictionary
)
stack = monet.plotMeanGenotypeStack(aggData, styleS)
monet.quickSaveFigure(
    stack,
    basePath + experimentString + "_Stack.png", format="png"
)
#############################################################################
landscapeData = monet.loadLandscapeData(filenames)
aggregatedNodesData = monet.aggregateGenotypesInLandscape(
    landscapeData,
    aggregationDictionary
)
geneSpatiotemporals = monet.getGenotypeArraysFromLandscape(aggregatedNodesData)
genes = geneSpatiotemporals["genotypes"]
plotsArray = monet.plotGenotypeArrayFromLandscape(
    geneSpatiotemporals,
    style={"aspect": 12, "cmap": cmaps}
)
for i in range(0, len(genes)):
    geneIndex = i
    fig = plotsArray["plots"][geneIndex]
    monet.quickSaveFigure(
        fig,
        basePath + experimentString + "_Heat" + genes[i] +  ".png"
    )

overlay = monet.plotGenotypeOverlayFromLandscape(
    geneSpatiotemporals,
    style={"aspect": 10, "cmap": cmaps}
)
monet.quickSaveFigure(overlay, basePath + experimentString + "HeatO.png")
