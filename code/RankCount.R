library(combinat)    ### Used for hcube
library(igraph)      ### Used for the adjacency graph
library(matrixStats) ### Used for logSumExp

### Constant scalars, vectors and matrices to be used in the computations
MAXSIZE = 63
MAXDIST = 800
MAXITER = 1e7
RANGE = 0:MAXSIZE
LOG_CATALAN = lchoose(2 * RANGE, RANGE) - log(RANGE + 1)
LOG_BINOM = outer(0:MAXDIST, 0:MAXDIST, lchoose)
LOG_CENTRAL = LOG_BINOM[cbind(RANGE + 2, floor((RANGE + 1)/2) + 1)]
LOG_CYCLENUMS = RANGE * log(RANGE + 2)

### Computes the multinomial coefficient [L1 + L2 + ... + LK ; L1, L2, ..., LK] with the numbers in the input list
### If offset is non-zero, adds it to the numerator; if logSpace = TRUE, the result is produced in log-space.
computeMultinomial = function(List, offset = 0, logSpace = TRUE) {
  partSums = cumsum(List)
  curSum = sum(LOG_BINOM[cbind(offset + partSums + 1, List + 1)])
  if (!logSpace) {
    curSum = exp(curSum)
  }
  curSum
}

### Computes the number of optimal rank scenarios between two genomes, using the recursion in Pereira & Meidanis
computeNumScenarios = function(maxK = 10, logSpace = TRUE) {
  op = ifelse(logSpace, magrittr::add, magrittr::multiply_by)
  combine = ifelse(logSpace, matrixStats::logSumExp, sum)
  remove = ifelse(logSpace, magrittr::subtract, magrittr::divide_by)
  transform = ifelse(logSpace, log, identity)
  minL = floor(maxK/2) + 1
  auxZ = 1:minL
  Z = remove(ifelse(rep(logSpace, minL), (auxZ - 2) * log(auxZ), auxZ^(auxZ - 2)), transform(factorial(auxZ - 1)))
  Mat = matrix(transform(0), maxK + 1, maxK + 1)
  Mat[1,1] = transform(1)
  for (k in 1:maxK) {
    for (l in 1:maxK) {
      if (l <= k && l > floor(k/2)) {
        iRange = 0:(ceiling(k/2)-1)
        jRange = 0:(l-1)
        term1 = op(Mat[2 * iRange + 1, jRange + 1, drop = FALSE], Mat[k - 2 * iRange, l - jRange, drop = FALSE])
        Mat[k + 1, l + 1] = combine(c(Mat[k + 1, l + 1], term1))
        maxZ = min(l, (ceiling(k/2) - 1))
        if (maxZ > 0) {
          zRange = 1:maxZ
          mults = transform(ceiling(k/2) - zRange)
          term2 = op(op(Mat[cbind(k - 2 * zRange + 1, l - zRange + 1)], mults), Z[zRange])
          Mat[k + 1, l + 1] = combine(c(Mat[k + 1, l + 1], term2))
        }
        Mat[k + 1, l + 1] = remove(Mat[k + 1, l + 1], transform(l))
      }
    }
  }
  for (l in 0:maxK) {
    Mat[, l + 1] = op(Mat[, l + 1], transform(factorial(l)))
  }
  Mat
}

signedPermToAdjGraph = function(signedGenome, originalChrLengths) {
  fullSignedGenome = unlist(signedGenome)
  n = length(fullSignedGenome)
  stopifnot(n >= 2)
  stopifnot(all(range(abs(fullSignedGenome)) == c(1, n)))
  tails = 1:n
  heads  = n + (1:n)
  allPairs = matrix(NA, 2 * n, 2)
  pos = 0
  val = 0
  for (index1 in 1:length(originalChrLengths)) {
    curLength = originalChrLengths[index1]
    if (curLength > 1) {
      firsts  = heads[val + (1:curLength)][-curLength]
      seconds = tails[val + (1:curLength)][-1]
      allPairs[pos + (1:(curLength - 1)), ] = cbind(firsts, seconds)
    }
    val = val + curLength
    pos = pos + curLength
  }
  for (index2 in 1:length(signedGenome)) {
    curChromosome = signedGenome[[index2]]
    curLength = length(curChromosome)
    if (curLength > 1) {
      firsts  = ifelse(curChromosome > 0, heads[abs(curChromosome)], tails[abs(curChromosome)])[-curLength]
      seconds = ifelse(curChromosome > 0, tails[abs(curChromosome)], heads[abs(curChromosome)])[-1]
      allPairs[pos + (1:(curLength - 1)), ] = cbind(firsts, seconds)
    }
    pos = pos + curLength
  }
  allPairs = allPairs[!is.na(allPairs[,1]) & !is.na(allPairs[,2]), ]
  G = graph(t(allPairs), directed = FALSE)
  G
}

adjGraphToCompLengths = function(G) {
  CL = clusters(G)
  allPathLengths = c()
  allCycleLengths = c()
  for (index in 1:CL[[3]]) {
    curSubgraph = induced_subgraph(G, which(CL[[1]] == index))
    curSize = ecount(curSubgraph)
    if (vcount(curSubgraph) == curSize) {
      allCycleLengths = c(allCycleLengths, curSize)
    } else {
      allPathLengths = c(allPathLengths, curSize)
    }
  }
  allCycleLengths = allCycleLengths[allCycleLengths != 2]
  allPathLengths  = allPathLengths[allPathLengths != 0]
  output = list(cycles = allCycleLengths, paths = allPathLengths)
  output
}

convertGenomePairsToLists = function(Dir = "eut", outputFile = "eut_k_lists.txt") {
  initDir = getwd()
  setwd(Dir)
  LF = list.files()
  outputLines = rep("", length(LF))
  for (index in 1:length(LF)) {
    inputFile = LF[index]
    print(paste("Currently processing", inputFile))
    f = file(inputFile, 'r')
    Lines = readLines(f)
    close(f)
    Lines = Lines[!grepl("^\\#", Lines)]
    Lines = Lines[Lines != ""]
    startInds = which(grepl("^>", Lines))
    startLines = Lines[startInds]
    startLines = sapply(strsplit(startLines, "\t"), function(x) {x[1]})
    startLines = gsub("^>", "", startLines)
    stopifnot(length(startLines) == 2)
    perm1 = Lines[(startInds[1] + 1):(startInds[2]-1)]
    perm1 = gsub(" \\$$", "", perm1)
    perm1 = lapply(perm1, function(x) { as.integer(unlist(strsplit(x, " "))) })
    originalChrLengths = sapply(perm1, length)
    L1 = sum(originalChrLengths)
    stopifnot(all(unlist(perm1) == 1:L1))
    perm2 = Lines[(startInds[2] + 1):length(Lines)]
    perm2 = gsub(" \\$$", "", perm2)
    perm2 = lapply(perm2, function(x) { as.integer(unlist(strsplit(x, " "))) })
    stopifnot(sum(sapply(perm2, length)) == L1)
    curGraph = signedPermToAdjGraph(perm2, originalChrLengths)
    res = adjGraphToCompLengths(curGraph)
    outputLines[index] = paste0(startLines[1], " ", startLines[2], " [", paste(res$cycles, collapse = ", "), "]", 
                                                                   " [", paste(res$paths,  collapse = ", "), "]")
  }
  setwd(initDir)
  f = file(outputFile, 'w')
  writeLines(outputLines, f)
  close(f)
}

### This function reads in the lists of cycle and path lengths for each genome pair
readLists = function(inputFile) {
  f = file(inputFile, 'r')
  Lines = readLines(f)
  close(f)
  L = length(Lines)
  List = c()
  for (curLine in Lines) {
    curLine = gsub("hg19 ", "", curLine)
    curInds = regexpr("[a-zA-Z0-9]*\\ ", curLine)
    curName = substr(curLine, curInds, curInds + attr(curInds, "match.length") - 2)
    curLine = gsub(paste(curName, ""), "", curLine)
    curLine = strsplit(curLine, "\\]\\ \\[")
    curLine = unlist(curLine)
    cycles = substr(curLine[1], 2, nchar(curLine[1]))
    paths = substr(curLine[2], 1, nchar(curLine[2]) - 1)
    cycles = gsub(" ", "", cycles)
    paths  = gsub(" ", "", paths)
    cycles = as.numeric(unlist(strsplit(cycles, ",")))
    paths  = as.numeric(unlist(strsplit(paths,  ",")))
    List[[curName]] = list(cycles = cycles, paths = paths)
  }
  List
}

### This function counts the log number of intermediates from the lists of cycle and path lengths in each genome pair
countIntermeds = function(List) {
 L = length(List)
 numIntermeds = rep(0, L)
 names(numIntermeds) = names(List)
 for (ind in 1:L) {
   curList = List[[ind]]
   curCycles = curList$cycles
   stopifnot(all(curCycles %% 2 == 0))
   curPaths  = curList$paths
   term1 = sum(LOG_CATALAN[(curCycles/2) + 1])
   term2 = sum(LOG_CENTRAL[curPaths + 1])
   numIntermeds[ind] = term1 + term2
 }
 numIntermeds
}

### This function counts the log number of optimal scenarios from a list of cycle and path lengths in genome pairs.
### Note that when the total number of combinations exceeds MAXITER, the last few dimensions are clamped to extreme
### values, meaning that the result is a lower bound and an upper bound (always in log space) which may be different.
countScenarios = function(List) {
  LOG_SCENARIOS = computeNumScenarios(maxK = MAXSIZE, logSpace = TRUE)
  L = length(List)
  upper = rep(0, L)
  lower = rep(0, L)
  names(upper) = names(List)
  names(lower) = names(List)
  for (index in 1:L) {
    curName = names(List)[index]
    print(curName)
    curList = List[[index]]
    curCycles = sort(curList$cycles)
    stopifnot(all(curCycles %% 2 == 0))
    halfCycles = curCycles/2
    shiftedHalfCycles = halfCycles - 1
    cycleSum = sum(shiftedHalfCycles)
    stopifnot(shiftedHalfCycles[1] >= 1)
    startTerm = sum(LOG_CYCLENUMS[shiftedHalfCycles])
    curPaths  = sort(curList$paths, decreasing = TRUE)
    numPaths = length(curPaths)
    stopifnot(numPaths >= 2)
    range1 = (floor(curPaths/2)[1] + 1):curPaths[1]
    N1 = length(range1)
    range2 = (floor(curPaths/2)[2] + 1):curPaths[2]
    N2 = length(range2)
    outR = as.vector(outer(range1, range2, "+"))
    repR2 = rep(range2, each = N1)
    sceMat = outer(LOG_SCENARIOS[cbind(curPaths[1], range1) + 1], LOG_SCENARIOS[cbind(curPaths[2], range2) + 1], "+")
    if (numPaths < 3) {
      startCoeff = computeMultinomial(shiftedHalfCycles)
      coeffMat = LOG_BINOM[cbind(cycleSum + range1, range1) + 1] + LOG_BINOM[cbind(cycleSum + outR, repR2) + 1]
      sum = startTerm + startCoeff + matrixStats::logSumExp(coeffMat + sceMat)
      curLower = sum
      curUpper = sum
      lower[index] = curLower
      upper[index] = curUpper
      next
    }
    curDimensions = ceiling(curPaths/2)
    print(curDimensions)
    preRange = 1:2 ### these arguments are matricized
    curProds = cumprod(curDimensions[-preRange])
    curFirst = 3   
    curLast = max(which(curProds < MAXITER)) + 2 
    curRange = curFirst:curLast ### these arguments are fully explored
    postRange = c()
    if (numPaths > curLast) {
      postRange = (curLast + 1):numPaths ### these arguments are clamped
    }
    exact = (length(postRange) == 0)
    if (exact) {
      print("Amenable to exact computation!")
    } else {
      print(paste("Setting the last", length(postRange), "dimensions to their extreme values"))
    }
    partProd = curProds[length(curRange)]
    curStarts = floor(curPaths/2)
    curLower = -Inf
    curUpper = -Inf
    upperVals = c()
    lowerVals = c()
    if (!exact) {
      upperVals = curPaths[postRange]
      lowerVals = floor(upperVals/2) + 1
    }
    upperSum = sum(upperVals)
    lowerSum = sum(lowerVals)
    myCube = combinat::hcube(curDimensions[curRange])
    print(paste("There are", partProd, "rows to process"))
    overallMult = ifelse(exact, 0, sum(log(curDimensions[postRange])))
    for (ind in 1:nrow(myCube)) {
      if (ind %% 1e5 == 0) {
        print(paste("Processing row", ind))
      }
      curRow = myCube[ind,] + curStarts[curRange]
      curProd = sum(LOG_SCENARIOS[cbind(curPaths[curRange] + 1, curRow + 1)])
      curProdUpper = curProd + ifelse(exact, 0, sum(LOG_SCENARIOS[cbind(curPaths[postRange], upperVals) + 1]))
      curProdLower = curProd + ifelse(exact, 0, sum(LOG_SCENARIOS[cbind(curPaths[postRange], lowerVals) + 1]))
      curSum = sum(curRow) + cycleSum
      curCoeff = computeMultinomial(c(curRow, shiftedHalfCycles))
      curCoeffUpper = curCoeff + ifelse(exact, 0, computeMultinomial(upperVals, offset = curSum))
      upperSumF = curSum + upperSum
      coeffMatUpper = LOG_BINOM[cbind(upperSumF + range1, range1) + 1] + LOG_BINOM[cbind(upperSumF + outR, repR2) + 1]
      curCoeffLower = curCoeff + ifelse(exact, 0, computeMultinomial(lowerVals, offset = curSum))
      lowerSumF = curSum + lowerSum
      coeffMatLower = LOG_BINOM[cbind(lowerSumF + range1, range1) + 1] + LOG_BINOM[cbind(lowerSumF + outR, repR2) + 1]
      curTermUpper = curProdUpper + curCoeffUpper + matrixStats::logSumExp(coeffMatUpper + sceMat)
      curTermLower = curProdLower + curCoeffLower + matrixStats::logSumExp(coeffMatLower + sceMat)
      curUpper = matrixStats::logSumExp(c(curUpper, curTermUpper))
      curLower = matrixStats::logSumExp(c(curLower, curTermLower))
    }
    curLower = curLower + startTerm + overallMult
    curUpper = curUpper + startTerm + overallMult
    lower[index] = curLower
    upper[index] = curUpper
  }
  output = rbind(lower, upper)
  output
}

### This function reads in a specified input and produces the distances and number of intermediates and scenarios
mainDriver = function(inputDir = "eut", outFile = "eut_k_lists.txt") {
  result = convertGenomePairsToLists(Dir = inputDir, outputFile = outFile)
  Lists = readLists(outFile)
  dists = sapply(Lists, function(x) {sum(x$cycles) - 2 * length(x$cycles) + sum(x$paths)})
  numIntermeds = exp(countIntermeds(Lists))
  numScenarios = countScenarios(Lists)
  exactScenarios = which(abs(numScenarios[1,] - numScenarios[2,]) < log(10))
  numScenarios[, exactScenarios]  = signif(exp(numScenarios[, exactScenarios]), 3)
  if (length(exactScenarios) < length(Lists)) {
    numScenarios[, -exactScenarios] = round(numScenarios[, -exactScenarios] / log(10))
  }
  print(paste("The numbers of scenarios for", paste(names(exactScenarios), collapse = ", "), "are exact"))
  print("For the remaining ones, the values is the exponent of the closest power of 10")
  output = cbind(d = dists, scenariosLB = numScenarios[1,], scenariosUB = numScenarios[2,], inter = numIntermeds)
  output
}
