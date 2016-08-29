// If you want to look up multiple patients at once and create *one* output file, use the lookUpMultiplePatients() function instead of lookUpOnePatient()
// This does not support adding in a 'vars' file and is generally not too useful as you'll probably want seperate output files for seperate patients
// E.g. use: 
// 		lookUpMultiplePatients([['cardiomyopathy', ['HP:0001637']], ['Orthostatic hypotension', ['HP:0001278']]])

if (process.argv.length < 6 || (process.argv[5] != 'x' && process.argv.length < 11)) {
	console.log('\n******************************************************************************************************************')
	console.log('\nPlease provide the following arguments:')
	console.log('$ node prioritization.js path predictions db vars [type] [cutoff] [results] [causalgene] [ID] HPO1 HPO2 HPO3 etc.')
	console.log('\n\tpath: the directory in which to store the output files (end with /)')
	console.log("\tpredictions: the name of your output file containing all prioritized genes (don't end with .txt)")
	console.log('\tdb: the location of your database, e.g. /Volumes/Promise_RAID/GeneNetwork/Sipko/05-2016/Tessa/hpodb_11jul')
	console.log("\n\tvars: file containing relevant variants (after filtering for MAF/variant impact/whatever else you fancy)")
	console.log("\t      (don't want to match with relevant variants? use 'x' instead)")
	console.log("\ttype: 'vcf' if vars is a VEP annotated vcf file taken through Patrick's script, \n\t      'cg2' (duo) or 'cg3' (trio) if cartagenia output from Kristin")
	console.log("\t      'gav' if vars is Gavin-output")
	console.log("\tcutoff: the number of top prioritized genes to overlap the vcf file with")
	console.log("\t      (only provide cutoff if vars is not 'x'!)")
	console.log("\tresults: the name of the output file in which to store the gene shortlist (don't end with .txt)")
	console.log("\t      (only provide results if vars is not 'x'!)")
	console.log("\tcausalgene: if the causal gene is known, enter here (in the 'ENSG00000158055' format) to report its shortlist rank\n\t      (only provide causalgene if vars is not 'x', no causal gene known? put 'x' here as well")
	console.log("\tID: add a number or other identifier (or 'x' again) here, very useful if you're running many of these sequentially\n\t      (only provide ID if vars is not 'x'!)")
	console.log("\nDO NOT FORGET: \nThen add all relevant HP terms (space-delimited) in the 'HP:0000482' format.")
	console.log('\nLooking up multiple patients (for one output file) is an option too, have a look at the comments at the top of this script\n')
	console.log('\n***********************************************************************************************************************')

} else {

	if (process.argv.length > 10 ) {
		console.log('ID: ' + process.argv[10])
	}

	var async = require('async')
	var fs = require('fs')
	var _ = require('lodash')
	var level = require('level')

	var termDBLocation = process.argv[4]

	var termDB = level(termDBLocation, {valueEncoding: 'binary'})

	var geneData = fs.readFileSync('files/ENSGToGeneNameHGNCBiotypeChromosomeStartStopStrandAndDescriptionV75.txt.filtered.txt', 'utf8')

	var filename = process.argv[3]
	var path = process.argv[2]
	var outputFile = path + filename + '.txt'
	var wstream = fs.createWriteStream(outputFile)

	if (process.argv[5] != 'x') {	// is a variants file included?
		var type = process.argv[6]
		var hpoArgs = process.argv.slice(11)
		var cutoff = process.argv[7]
		var resultsFile = path + process.argv[8] + '_' + cutoff + '.txt'
		var streamMatch = fs.createWriteStream(resultsFile)
		var vars = fs.readFileSync(process.argv[5], 'utf8')
	} else {
		var hpoArgs = process.argv.slice(6)
	}

	lookUpMultiplePatients = function (data) {
		wstream.write('Based on the HPO terms: ' + data + '\n')
		wstream.write('Patient\tPrioritized genes (gene name|gene code|z score|gene rank)\n')

		var totalPatientsInfo = data

		counter = 0

		async.mapSeries(totalPatientsInfo, function (item, callback) {		// for each patient


			counter += 1
			var patient = item[0]
			var HPOs = item[1]

			console.log('processing patient ' + counter + ' of ' + totalPatientsInfo.length)

			combineTerms(HPOs, function(err, results) {			
			// results contains an array of objects ({gene: geneN, code: geneC, score: z, }), sorted according to the score

				wstream.write(patient + '\t')

				for (i = 0; i < results.length; i++) {			// for each gene in HPO zScore list
					resultLine = results[i]
					var rank = i + 1				
					wstream.write(resultLine['gene'] + '|' + resultLine['code'] + '|' + resultLine['score'] + '|' + rank + '\t' )
				}

				wstream.write('\n') // end patient row in file

				callback(err, results) 
			})
		})
	}

	lookUpOnePatient = function (data) {
		wstream.write('Based on the HPO terms: ' + data + '\n')
		wstream.write('Gene name\tGene code\tZ score\tRank\n')

		var HPOs = data

		if (process.argv[5] != 'x') {
			streamMatch.write('Based on the HPO terms: ' + data + '\n')
			if (type == 'vcf') {
				streamMatch.write('Shortlist Rank\tGene name\tGene code\tZ score\tPrediction rank\tImpact\n')			
			} else if (type == 'cg2' || type == 'cg3') {
				streamMatch.write('Shortlist Rank\tGene name\tGene code\tZ score\tPrediction rank\n')
				// streamMatch.write('Shortlist Rank\tGene name\tGene code\tZ score\tPrediction rank\tCoding effect\n')
			} else if (type == 'gav') {
				streamMatch.write('Shortlist Rank\tGene name\tGene code\tZ score\tPrediction rank\n')
			}
			
			var varLines = vars.split('\n')
			var varArray = []

			if (type == 'gav') {
				for (i = 0; i < varLines.length; i++) {
					var varLine = varLines[i].split('\t')
					// console.log(varLine)
					// console.log(varLine[7])
					if (varLine[7] !== undefined) {
						var varInfo = varLine[7].split(';')
						// console.log(varInfo)
						for (j = 0; j < varInfo.length; j++) {
							// console.log(varInfo[j].substring(0,3))
							if (varInfo[j].substring(0,3) == 'RLV') {
								// console.log(varInfo[j])
								var varRLV = varInfo[j].split('|')
								// console.log(varRLV)
								varArray.push(varRLV)
							}
						}	
					}	
				} //console.log(varArray)

			} else {
				for (i = 0; i < varLines.length; i++) {
					var varLine = varLines[i].split('\t')
					varArray.push(varLine)
				}
			}
		}

		combineTerms(HPOs, function(err, results) {	

		// results contains an array of objects ({gene: geneN, code: geneC, score: z, }), sorted according to the score

			var previousGene = undefined
			var k = 1

			for (i = 0; i < results.length; i++) {			// for each gene in HPO zScore list
				resultLine = results[i]
				var rank = i + 1
				wstream.write(resultLine['gene'] + '\t' + resultLine['code'] + '\t' + resultLine['score'] + '\t' + rank + '\n' )
				if (resultLine['code'] == process.argv[9]) {
					console.log('Causal gene:', resultLine['gene'], '-', resultLine['code'])
					console.log('Predicted at rank ' + rank)
				}

				if (process.argv[5] != 'x') {
					if (i < cutoff) {
						// console.log(i)
						for (j = 0; j < varArray.length; j++) {
							// console.log('.')
							if (type == 'vcf') {
								if (varArray[j][1] == resultLine['code']) {
									if (previousGene != resultLine['code']) {
										streamMatch.write(k + '\t' + resultLine['gene'] + '\t' + resultLine['code'] + '\t' + resultLine['score'] + '\t' + rank + varArray[j][3] + '\n')
										previousGene = resultLine['code']
										if (process.argv[9] != 'x') {
											if (resultLine['code'] == process.argv[9]) {
												// console.log('Patient: ' + varArray[j][0])
												// console.log('Causal gene:', resultLine['gene'], '-', resultLine['code'])
												console.log('Shortlist rank ' + k)	
											}
										}
										k += 1								
									}
								}
							} else if (type == 'cg2') {
								if (varArray[j][107] == resultLine['code'] && varArray[j].length != 1) {
									if (previousGene != resultLine['code']) {
										streamMatch.write(k + '\t' + resultLine['gene'] + '\t' + resultLine['code'] + '\t' + resultLine['score'] + '\t' + rank  + '\n')
										previousGene = resultLine['code']
										if (process.argv[9] != 'x') {
											if (resultLine['code'] == process.argv[9]) {
												// console.log('Causal gene:', resultLine['gene'], '-', resultLine['code'])
												console.log('Shortlist rank ' + k)	
											}
										}
										k += 1
									}
								}

							} else if (type == 'cg3') {
								//console.log(varArray[j][126])
								if (varArray[j][126] == resultLine['code'] && varArray[j].length != 1) {
									if (previousGene != resultLine['code']) {
										streamMatch.write(k + '\t' + resultLine['gene'] + '\t' + resultLine['code'] + '\t' + resultLine['score'] + '\t' + rank + '\n')
										previousGene = resultLine['code']
										if (process.argv[9] != 'x') {
											if (resultLine['code'] == process.argv[9]) {
												// console.log('Causal gene:', resultLine['gene'], '-', resultLine['code'])
												console.log('Shortlist rank ' + k)	
											}
										}
										k += 1
									}
								}

							} else if (type == 'gav') {
								// console.log('type:gav')
								// console.log(varArray[j])
								// console.log(varArray[j][3], resultLine['gene'])
								if (varArray[j][2] == resultLine['gene']) {
									//console.log('yay, match')
									if (previousGene != resultLine['gene']) {
										// var b = varArray[j].length - 4
										// streamMatch.write(k + '\t' + resultLine['gene'] + '\t' + resultLine['code'] + '\t' + resultLine['score'] + '\t' + rank + '\t' + varArray[j][b] + '\n')
										streamMatch.write(k + '\t' + resultLine['gene'] + '\t' + resultLine['code'] + '\t' + resultLine['score'] + '\t' + rank + '\n')
										previousGene = resultLine['gene']
										if (process.argv[9] != 'x') {
											if (resultLine['code'] == process.argv[9]) {
												// console.log('Causal gene:', resultLine['gene'], '-', resultLine['code'])
												console.log('Shortlist rank ' + k)	
											}
										}
										k += 1								
									}
								}

							} else {
								console.log('Vars type not specified correctly!')
							}
							// if ( indexOf(resultLine['code']) Object.keys(varObject) )
						}
					}
				}
			}
			if (process.argv[5] != 'x') {
				var kEnd = k -1
				console.log('Total length of shortlist: ' + kEnd + '.\n')
			}
		})
	}


	combineTerms = function (HPOArray, callback) {
	// called for each patient (+ set of phenotypes): all genes are prioritzed
	// then looks up the rank for each involved gene in HPOGeneArray

		var resultsArray = []
		async.mapSeries(HPOArray, function(hpoTerm, cb) {		// for each HPO term

			dbLookUp(hpoTerm, function(err, dbResult) {		// look up all genes & z-scores for that term (= dbResult)
				return cb(null, dbResult)
			})
		}, function(err, results) {
			// now results is an array of {gene: geneN, code: geneC, score: z, } for each gene for that HPO term
			results = _.compact(results)
			var resultsArray = []

			for (i = 0; i < results.length; i++) {
				if (resultsArray.length === 0) {							// first HPO term: resultsArray is result dbLookUp
					resultsArray = results[i]
				} else {								// next terms: add result dbLookup to scores in resultsArray (sum results for different phenotypes)
					for (j = 0; j < resultsArray.length; j++) {
						resultsArray[j]['score'] += results[i][j]['score']
					}
				}
			}
			// now the zScores are summed in resultsArray for all HPO terms

			resultsArray = resultsArray.sort(function (a, b) {	// sort array of summed z-scores & genes
				return b.score - a.score
			})
			callback(err, resultsArray)
		})
	} 



	dbLookUp = function (HPOterm, callback) {
	// Looks up and returns gene prioritization for one HPO term (all genes and their Z-scores, unsorted)
		
		var genesList = []
		var geneNames = []
		lines = geneData.split('\n')

		for (i = 1; i < lines.length; i++) {
			if (lines[i] !== "") {
				var line = lines[i]
				var splitLine = line.split('\t')
				var geneCode = splitLine[0]	// ENSG code of the gene
				var geneName = splitLine[1] // name of the gene
				genesList.push(geneCode)
				geneNames.push(geneName)
			}
		}

		var key = 'RNASEQ!PREDICTIONS!HPO!' + HPOterm
		termDB.get(key, function(err, value) {		// value is a buffer of 54,000 Z-scores

			if (err) {
				if (err.name == "NotFoundError") {
					console.log("not found: " + HPOterm)
					callback(null, null)
				} else {
					return callback(err)
				}

			} else {		// no error
				var geneZScoresArray = []
				for (var i = 1; i < value.length / 2; i++) { // for all z-scores in db (one per gene)
					var geneC = genesList[i - 1]
					var geneN = geneNames[i - 1]
					if (geneN != undefined) {
						var z = (value.readUInt16BE(i * 2) - 32768) / 1000
						geneZScoresArray.push({gene: geneN, code: geneC, score: z, })
					}
				}
				callback(null, geneZScoresArray)
			}
		})
	}

	lookUpOnePatient(hpoArgs)
	//E.g.:
		
		// lookUpOnePatient(['HP:0001637']) // cardiomyopathy
		// lookUpOnePatient(['HP:0001278']) // Orthostatic hypotension
	//	lookUpOnePatient(["HP:0002987","HP:0010628","HP:0003676","HP:0003560","HP:0000007","HP:0003691"]) // muscular dystrophy (LG)


	// lookUpOnePatient(["HP:0000505","HP:0000482","HP:0000483","HP:0007663","HP:0000006","HP:0007676"])

	// lookUpOnePatient(["HP:0001874", "HP:0001419", "HP:0002718", "HP:0004313", "HP:0000951"]) // WAS
}