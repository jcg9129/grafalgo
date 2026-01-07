/** @file setCoverRandom.mjs
 *
 *  @author Jon Turner
 *  @date 2024
 *  This is open source software licensed under the Apache 2.0 license.
 *  See http://www.apache.org/licenses/LICENSE-2.0 for details.
 */

import {assert, EnableAssert as ea} from '../../common/Assert.mjs';
import List from '../../dataStructures/basic/List.mjs';
import Graph from '../../dataStructures/graphs/Graph.mjs';
import { shuffle } from '../../common/Random.mjs';
import { randomBigraph, regularize }
		from '../../graphAlgorithms/misc/RandomGraph.mjs';
import { randomFraction } from '../../common/Random.mjs';
import setCoverSizeBound from './setCoverSizeBound.mjs';
import setCoverSplitBound from './setCoverSplitBound.mjs';
import setCoverLabelBound from './setCoverLabelBound.mjs';

/** Generate a random set cover instance.
 *  @param k is the number of sets
 *  @param h is the number of elements in the base set
 *  @param coverage is the number of times each item appears in a subset
 *  @param randomWeight is a function that returns a random number
 *  @param args collects remaining arguments into an array of arguments
 *  for randomWeight
 *  @return tuple [g,weight,lowerBounds,upperBound] and lowerBounds is an
 *  array of three lower bounds on the opimimum cover weight and upperBound
 *  is an upperBound on the weight, specifically it is the weight of
 *  a secret cover that is embedded in the constructed instance.
 */
export default function setCoverRandom(k, h, coverage, randomWeight, ...args) {
	let subSize = coverage*h/k; // average subset size

	// allow for small variation in coverage
	let secretCoverage = (coverage <= 2.1 ? 1 : 1.05);
	let camoCoverage = coverage - secretCoverage;

	// determine number of subsets in secret and camouflage
	let secretWidth = Math.round(h*secretCoverage/subSize);
	let camoWidth = k - secretWidth

	// allow some irregularity in camo coverage 
	let camoRegularity = Math.max(1, Math.log2(camoCoverage));

	// create graphs for  secret and comouflage
	let secret = randomBigraph(secretWidth, h*secretCoverage/secretWidth, h);
	let items = new List(h); items.range(secretWidth+1,secretWidth+h);
		regularize(secret, secretCoverage, items); 
	items.range(1, secretWidth);
		regularize(secret, h*secretCoverage/secretWidth, items,
				   Math.max(1, Math.log2(subSize)));
	let camo = randomBigraph(camoWidth, h*camoCoverage/camoWidth, h);
	items.range(camoWidth+1,camoWidth+h);
		regularize(camo, camoCoverage, items, camoRegularity);
	items.range(1, camoWidth);
		regularize(camo, h*camoCoverage/camoWidth, items,
				   Math.max(1, Math.log2(subSize)));

	// combine the graphs
	let g = new Graph(k+h, subSize*k); g.setBipartition(k);
	for (let e = secret.first(); e; e = secret.next(e))
		g.join(secret.left(e), secret.right(e)+(k-secretWidth));
	for (let e = camo.first(); e; e = camo.next(e))
		g.join(camo.left(e)+secretWidth, camo.right(e)+secretWidth);

	// assign costs to subsets, with lower costs for subsets in secret
	let weightList = new Float32Array(k);
	for (let j in weightList) weightList[j] = randomWeight(...args);
	weightList.sort();
	let weight = new Float32Array(k+1); let upperBound = 0;
	let s = 1; let sr = secretWidth;	// next secret subset, # remaining
	let t = secretWidth+1; let tr = camoWidth; // next camo subset, # remaining
	for (let w of weightList) {
		if (sr && randomFraction() < Math.sqrt(coverage)*(sr/(sr+tr))) {
			weight[s++] = w; sr--; upperBound += w;
		} else {
			weight[t++] = w; tr--;
		}
	}
	// upperBound is now the total weight of the secret cover

	// now, scramble graph, while keeping outputs fixed
	let outputs = new Set();
	for (let i = k+1; i <= k+h; i++) outputs.add(i);
	let [vp] = g.scramble(outputs);
	shuffle(weight, vp.slice(0,k+1)); 

	g.sortAllEplists()

	let lowerBounds = [ setCoverSizeBound(g,weight),
						setCoverSplitBound(g,weight),
						setCoverLabelBound(g,weight) ];

	return [g, weight, lowerBounds, upperBound];
}
