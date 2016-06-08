package net.maizegenetics.gwas.imputation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.regex.Pattern;

public class TransitionProbabilityWithVariableRecombination extends TransitionProbability{
	double averageRecombinationRate;
	int[] ratePosition = null;
	double[] rateAtPosition = null;
	int numberOfPositions = 0;
	
	public TransitionProbabilityWithVariableRecombination(String chromosome) {
		readRateFile(chromosome);
	}

	@Override
	public void setNode(int node) {
		if (node <= 0) return;
		double rrr = getRelativeRecombinationRate(node);
		int n = probabilityOfATransition.length;
		adjustedProbability = new double[n][n];
		int segmentLength = positions[node] - positions[node - 1];
		double m;
		for (int row = 0; row < n; row++) {
			double offdiagsum = 0;
			for (int col = 0; col < n; col++) {
				if (col != row) {
					m = -Math.log(1 - 2 * probabilityOfATransition[row][col] * rrr) * segmentLength / avgSegmentLength / 2;
					adjustedProbability[row][col] = (1 - Math.exp(-2*m)) / 2;
					offdiagsum += adjustedProbability[row][col];
				}
			}
			adjustedProbability[row][row] = 1 - offdiagsum;
		}
	}
	
	private double getRelativeRecombinationRate(int node) {
		int start = Arrays.binarySearch(ratePosition, positions[node - 1]);
		int end = Arrays.binarySearch(ratePosition, positions[node]);
		if (start < 0) start = -start - 2;
		if (start < 0) start = 0;
		if (end < 0) end = -end - 1;
		if (end >= numberOfPositions) end = numberOfPositions - 1;
		double sum = 0;
		for (int i = start; i<=end; i++) sum += rateAtPosition[i];
		
		return sum/(end - start + 1) / averageRecombinationRate;
	}
	
	private boolean readRateFile(String chromosome) {
		String filename = "/Volumes/Macintosh HD 2/results/recombination study/nam/final.Panzea.consolidated.B/with.missing/recombination.rate.chr" + chromosome + ".txt";
		Pattern tab = Pattern.compile("\t");
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			int nlines = 0;
			while (br.readLine() != null) nlines++;
			br.close();
			int npts = nlines - 1; 
			double sumRates = 0;
			ratePosition = new int[npts];
			rateAtPosition = new double[npts];
			br = new BufferedReader(new FileReader(filename));
			br.readLine();
			for (int i = 0; i < npts; i++) {
				String[] data = tab.split(br.readLine());
				ratePosition[i] = (int) Double.parseDouble(data[0]);
				rateAtPosition[i] = Double.parseDouble(data[1]);
				sumRates += rateAtPosition[i];
			}
			br.close();
			averageRecombinationRate = sumRates/npts;
			numberOfPositions = npts;
			return true;
		} catch(IOException e) {
			e.printStackTrace();
			return false;
		}
	}
}
