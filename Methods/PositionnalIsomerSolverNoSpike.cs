/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using PeptidAce.Utilities;
using PeptidAce.Utilities.Methods;
using PeptidAce.Iso.Structures;
using PeptidAce.Utilities.Interfaces;

namespace PeptidAce.Iso.Methods
{
    /// <summary>
    /// DEVELOPMENT
    /// This class is a simple copy of PositionnalIsomerSolver dedicated into trying to deconvolute mixed precursor elution
    /// curves using only mixed raws (no spiked synthetic peptide run)
    /// 
    /// I was not particularly successfull so far. Its based on amino acid coverage of each peptide variant. In theory, it should 
    /// be possible. I suspect having peptide isomers of different sequences would help. Positionnal Isomers are an added difficulty here.
    /// </summary>
    public class PositionnalIsomerSolverNoSpike
    {
        public int nbMinFragments = 5;
        public int nbMaxFragments = 5;
        public double precTolPpm = 8;
        public double prodTolPpm = 20;
        public long precision     = 100;
        private DBOptions dbOptions;
        
        private Samples MixedSamples;
        private Result mixedResult;
        public Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> characterizedPeptides;
        public Dictionary<Sample, List<MixedPrecursor>>         mixedPrecursors;

        public string OutputFolder { get { return dbOptions.OutputFolder + "CombinedNoSpike" + System.IO.Path.DirectorySeparatorChar;  } }
        private DBOptions CreateOptions(string fastaFile, string outputFolder, IConSol consol)
        {
            DBOptions dbOptions = new DBOptions(fastaFile, consol);
            dbOptions.precursorMassTolerance = new MassTolerance(precTolPpm, MassToleranceUnits.ppm);
            dbOptions.productMassTolerance = new MassTolerance(prodTolPpm, MassToleranceUnits.ppm);
            dbOptions.MaximumPeptideMass = 200000;
            dbOptions.OutputFolder = outputFolder;

            ProteaseDictionary proteases = ProteaseDictionary.Instance;
            //dbOptions.DigestionEnzyme = proteases["no enzyme"];
            dbOptions.DigestionEnzyme = proteases["top-down"];
            dbOptions.NoEnzymeSearch = false;
            dbOptions.DecoyFusion = false;
            dbOptions.MaximumNumberOfFragmentsPerSpectrum = 400;
            dbOptions.ToleratedMissedCleavages = 200;
            dbOptions.MinimumPeptideLength = 5;
            dbOptions.MaximumPeptideLength = 300;

            List<Modification> fixMods = new List<Modification>();
            dbOptions.fixedModifications = fixMods;

            List<Modification> varMods = new List<Modification>();
            foreach (string strMod in ModificationDictionary.Instance.Keys)
                varMods.Add(ModificationDictionary.Instance[strMod]);

            dbOptions.maximumVariableModificationIsoforms = 1024;
            dbOptions.variableModifications = varMods;

            dbOptions.addFragmentLoss = false;
            dbOptions.addFragmentMods = false;//Gives very bad results... might by off
            
            dbOptions.SaveMS1Peaks = true;
            dbOptions.SaveMSMSPeaks = true;
            dbOptions.LoadSpectraIfFound = true;

            dbOptions.NbPSMToKeep = 100;
            return dbOptions;
        }

        public void Solve(string[] mixedRaws, string fastaFile, string folderToOutputTo, IConSol conSol)
        {
            dbOptions = CreateOptions(fastaFile, folderToOutputTo, conSol);

            MixedSamples = new Samples(dbOptions);
            for (int i = 0; i < mixedRaws.Length; i++)
                MixedSamples.Add(new Sample(i + 1, 1, 1, mixedRaws[i], mixedRaws[i], 0, ""));

            //Precompute Mixed peptide identifications
            mixedResult = Ace.Start(dbOptions, MixedSamples, false, false);

            conSol.WriteLine("Computing gradient descents...");

            //Compute all usable spiked peptides
            characterizedPeptides = CharacterizedPrecursor.GetSpikedPrecursors(MixedSamples, mixedResult, dbOptions, nbMinFragments, nbMaxFragments);
            ExportSpikedSampleResult(characterizedPeptides, dbOptions);

            vsCSVWriter writerCumul = new vsCSVWriter(OutputFolder + "Results.csv");
            string titleCombined = "Mixed Sample,Precursor";            
            string curveStr = "Polynomial Curve,";
            string spikedIntensityStr = "Area under the curve,";
            foreach(double precursor in characterizedPeptides.Keys)
                foreach (CharacterizedPrecursor charPrec in characterizedPeptides[precursor].Values)
                {
                    titleCombined += "," + charPrec.Peptide.Sequence + " Charge " + charPrec.Charge;

                    if (charPrec.eCurveIntensityCount.Coefficients != null && charPrec.eCurveIntensityCount.Coefficients.Length == 3)
                        curveStr += "," + charPrec.eCurveIntensityCount.Coefficients[0] + "x^2 + " + charPrec.eCurveIntensityCount.Coefficients[1] + "x" + charPrec.eCurveIntensityCount.Coefficients[2];
                    else
                        curveStr += ",NA";

                    spikedIntensityStr += "," + charPrec.eCurveIntensityCount.Area;
                }
            writerCumul.AddLine(titleCombined);
            writerCumul.AddLine(curveStr);
            writerCumul.AddLine(spikedIntensityStr);

            //mixedPrecursors = new Dictionary<Sample, Dictionary<double, MixedPrecursor>>();
            mixedPrecursors = new Dictionary<Sample, List<MixedPrecursor>>();

            foreach (Sample mixedSample in MixedSamples) 
                mixedPrecursors.Add(mixedSample, MixedPrecursor.GetMixedPrecursors(mixedSample, mixedResult, dbOptions, characterizedPeptides));

                //Get the list of precursors to characterize
            foreach (Sample mixedSample in MixedSamples)
            {
                foreach (double keyMz in characterizedPeptides.Keys)
                {
                    List<Dictionary<Peptide, MaxFlowElutionCurve>> listOfRatios = new List<Dictionary<Peptide, MaxFlowElutionCurve>>();
                    foreach(MixedPrecursor mPrec in mixedPrecursors[mixedSample])
                        if(mPrec.MZ == keyMz)
                        {
                            // Compute Max Flow for this precursor
                            Dictionary<Peptide, MaxFlowElutionCurve> ratios = GetRatiosNoSpikes(mPrec, precision);
                            listOfRatios.Add(ratios);

                            ExportMixedSampleResult(ratios, mixedSample, mPrec, keyMz, dbOptions);
                        }
                    /*
                    string resultStr = vsCSV.GetFileName(mixedSample.sSDF) + "," + keyMz;
                    foreach (double precursor in characterizedPeptides.Keys)
                    {
                        foreach (Peptide charPrec in characterizedPeptides[precursor].Values)
                        {
                            double cumulArea = 0.0;
                            foreach (Dictionary<Peptide, ElutionCurve> ratios in listOfRatios)
                                if (ratios.ContainsKey(charPrec))
                                    cumulArea += ratios[charPrec].Area;
                            resultStr += "," + cumulArea;
                        }
                    }
                    writerCumul.AddLine(resultStr);//*/
                }
            }
            writerCumul.WriteToFile();
        }

        private static void ExportMixedSampleResult(Dictionary<Peptide, MaxFlowElutionCurve> ratios, Sample mixedSample, MixedPrecursor mixedPrecursor, double keyMz, DBOptions dbOptions)
        {
            // Export results in a file
            vsCSVWriter writerRatio = new vsCSVWriter(dbOptions.OutputFolder + @"IndividualNoSpike\" + vsCSV.GetFileName_NoExtension(mixedSample.sSDF) + "_" + keyMz + "MZ_" + mixedPrecursor.Queries[0].spectrum.RetentionTimeInMin + "min.csv");
            string titleIndividual = "Scan time,Total Area";
            foreach (Peptide charPep in ratios.Keys)
                titleIndividual += "," + charPep.Sequence;
            writerRatio.AddLine(titleIndividual);

            string line = "Total," + mixedPrecursor.eCurveIntensityCount.Area;
            foreach (Peptide charPep in ratios.Keys)
                line += "," + ratios[charPep].eCurvePerMs.Area;
            writerRatio.AddLine(line);

            for (int i = 0; i < mixedPrecursor.eCurveIntensityCount.intensityCount.Count; i++)
            {
                line = mixedPrecursor.eCurveIntensityCount.time[i] / (1000.0 * 60.0) + "," + mixedPrecursor.eCurveIntensityCount.intensityCount[i];
                foreach (Peptide charPep in ratios.Keys)
                    line += "," + ratios[charPep].eCurvePerMs.InterpolateIntensity(mixedPrecursor.eCurveIntensityCount.time[i]);
                writerRatio.AddLine(line);
            }
            writerRatio.WriteToFile();
        }

        private static void ExportSpikedSampleResult(Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> characterizedPeptides, DBOptions dbOptions)
        {
            foreach (double keyMz in characterizedPeptides.Keys)
            {
                foreach (Sample sample in characterizedPeptides[keyMz].Keys)
                {
                    vsCSVWriter writerRatio = new vsCSVWriter(dbOptions.OutputFolder + @"IndividualNoSpike\" + vsCSV.GetFileName_NoExtension(sample.sSDF) + "_" + keyMz + "MZ.csv");
                    string titleIndividual = "Scan time,Precursor Intensity,Intensity Per Millisecond";
                    foreach (ProductMatch pm in characterizedPeptides[keyMz][sample].AllFragments)
                        titleIndividual += "," + pm.Fragment.Name + pm.fragmentPos + "^" + pm.charge;
                    writerRatio.AddLine(titleIndividual);

                    foreach (Query query in characterizedPeptides[keyMz][sample].Queries)
                    {
                        string line = query.spectrum.RetentionTimeInMin + "," + query.spectrum.PrecursorIntensity + "," + query.spectrum.PrecursorIntensityPerMilliSecond;
                        foreach (ProductMatch pm in characterizedPeptides[keyMz][sample].AllFragments)
                        {
                            double intensity = 0.0;
                            foreach (ProductMatch pmSpec in query.psms[0].AllProductMatches)
                                if (pmSpec.charge == pm.charge && pmSpec.Fragment == pm.Fragment && pmSpec.fragmentPos == pm.fragmentPos)
                                    intensity = pmSpec.obsIntensity;
                            line += "," + intensity;
                        }
                        writerRatio.AddLine(line);
                    }
                    writerRatio.WriteToFile();
                }
            }
        }

        private Dictionary<Peptide, MaxFlowElutionCurve> GetRatiosNoSpikes(MixedPrecursor mixedPrecursor, long precision)        
        {
            Dictionary<Dictionary<Peptide, MaxFlowElutionCurve>, double> DicOfCurveErrors = new Dictionary<Dictionary<Peptide, MaxFlowElutionCurve>, double>();
            
            for (int nbProductsToKeep = nbMinFragments; nbProductsToKeep <= nbMaxFragments; nbProductsToKeep++)
            {
                bool validProducts = true;
                int nbIgnoredSpectrum = 0;
                if (validProducts)
                {
                    double cumulError = 0;
                    Dictionary<Peptide, MaxFlowElutionCurve> curves = new Dictionary<Peptide, MaxFlowElutionCurve>();
                    foreach (Query query in mixedPrecursor.Queries)
                    {
                        double timeInMilliSeconds = query.spectrum.RetentionTimeInMin * 60.0 * 1000.0;
  //                      double overFlow = 0;
                        double underFlow = 0;
                        double percentError = 0;
                        Dictionary<PeptideSpectrumMatch, SolvedResult> finalRatios = SolveFromFragmentScore(nbProductsToKeep, precision, query.spectrum.Peaks, dbOptions.productMassTolerance,
                                                                mixedPrecursor.eCurveIntensityPerMS.GetLocalArea(timeInMilliSeconds, timeInMilliSeconds + query.spectrum.InjectionTime),
                                                                query, out underFlow, out percentError, dbOptions.ConSole);

                        cumulError += underFlow;// percentError;
                        if (percentError < 0.5)
                        {
                            foreach (PeptideSpectrumMatch cPep in finalRatios.Keys)
                            {
                                if (!curves.ContainsKey(cPep.Peptide))
                                    curves.Add(cPep.Peptide, new MaxFlowElutionCurve(nbProductsToKeep));

                                //curves[cPep].AddPoint(timeInMilliSeconds, finalRatios[cPep].Ratio * query.spectrum.PrecursorIntensityPerMilliSecond);
                                curves[cPep.Peptide].eCurveCount.AddPoint(timeInMilliSeconds, finalRatios[cPep].Ratio * mixedPrecursor.eCurveIntensityCount.InterpolateIntensity(timeInMilliSeconds));
                                curves[cPep.Peptide].eCurvePerMs.AddPoint(timeInMilliSeconds, finalRatios[cPep].Ratio * mixedPrecursor.eCurveIntensityPerMS.InterpolateIntensity(timeInMilliSeconds));
                            }
                        }
                        else
                            nbIgnoredSpectrum++;

                        if (nbIgnoredSpectrum * 2 > mixedPrecursor.Queries.Count)
                            break;
                    }//End of foreach query

                    if (nbIgnoredSpectrum * 2 < mixedPrecursor.Queries.Count)
                    {
                        if (nbIgnoredSpectrum > 0)
                            Console.WriteLine("Ignored Spectrum : " + nbIgnoredSpectrum + " / " + mixedPrecursor.Queries.Count);

                        foreach (Peptide cPep in curves.Keys)
                            curves[cPep].Compute();

                        Dictionary<Peptide, MaxFlowElutionCurve> curvesToKeep = new Dictionary<Peptide, MaxFlowElutionCurve>();
                        foreach (Peptide cPep in curves.Keys)
                            if (curves[cPep].eCurvePerMs.Area > 0)
                                curvesToKeep.Add(cPep, curves[cPep]);

                        if (curvesToKeep.Count > 0)
                            DicOfCurveErrors.Add(curvesToKeep, cumulError);
                    }
                }
            }//End of for each nbProduct            

            Dictionary<Peptide, MaxFlowElutionCurve> averagedValues = ComputePeptideRatios(DicOfCurveErrors);
            return averagedValues;
        }

        public Dictionary<Peptide, MaxFlowElutionCurve> ComputePeptideRatios(Dictionary<Dictionary<Peptide, MaxFlowElutionCurve>,double> dicOfCurveErrorsP)
        {
            Dictionary<Dictionary<Peptide, MaxFlowElutionCurve>, double> dicOfCorrelations = new Dictionary<Dictionary<Peptide, MaxFlowElutionCurve>, double>();
            foreach (Dictionary<Peptide, MaxFlowElutionCurve> dicOfCurve in dicOfCurveErrorsP.Keys)
                dicOfCorrelations.Add(dicOfCurve, 1.0 / (double)dicOfCurveErrorsP.Count);

            //Purge worst curves
            Dictionary<Dictionary<Peptide, MaxFlowElutionCurve>, double> dicOfCurves = new Dictionary<Dictionary<Peptide, MaxFlowElutionCurve>, double>();
            if (dicOfCurveErrorsP.Count > 1)
            {
                double median = MathNet.Numerics.Statistics.Statistics.Median(dicOfCurveErrorsP.Values);
                double maxMed = median;// +0.5 * MathNet.Numerics.Statistics.Statistics.Variance(dicOfCurveErrorsP.Values);
                foreach (Dictionary<Peptide, MaxFlowElutionCurve> dic in dicOfCurveErrorsP.Keys)
                    if (dicOfCurveErrorsP[dic] <= maxMed)
                        dicOfCurves.Add(dic, dicOfCurveErrorsP[dic]);
            }
            else
                dicOfCurves = dicOfCurveErrorsP;

            Dictionary<Dictionary<Peptide, MaxFlowElutionCurve>, double> lastDicOfCurves = dicOfCurves;
            int nbRun = 2;
            while (nbRun > 0)
            {
                nbRun--;
                dicOfCurves = lastDicOfCurves;

                //Normalize already computed correlation factors for the remaning curves (sum must equal 1)
                double sumOfCorr = 0.0;
                foreach (Dictionary<Peptide, MaxFlowElutionCurve> dicOfCurve in dicOfCurves.Keys)
                    sumOfCorr += dicOfCorrelations[dicOfCurve];

                foreach (Dictionary<Peptide, MaxFlowElutionCurve> dicOfCurve in dicOfCurves.Keys)
                    dicOfCorrelations[dicOfCurve] /= sumOfCorr;

                //Compute average from weighted curves
                Dictionary<Peptide, double> average = new Dictionary<Peptide, double>();
                foreach (Dictionary<Peptide, MaxFlowElutionCurve> dicOfCurve in dicOfCurves.Keys)
                {
                    Dictionary<Peptide, double> areas = MixedPrecursor.GetAreas(dicOfCurve);
                    foreach (Peptide cPep in areas.Keys)
                    {
                        if (!average.ContainsKey(cPep))
                            average.Add(cPep, 0);
                        average[cPep] += areas[cPep] * dicOfCorrelations[dicOfCurve];
                    }
                }

                //Compute correlation between average and curves
                List<double> corrs = new List<double>();
                foreach (Dictionary<Peptide, MaxFlowElutionCurve> dicOfCurve in dicOfCurves.Keys)
                {
                    Dictionary<Peptide, double> elution = new Dictionary<Peptide, double>();
                    foreach (Peptide cPep in average.Keys)
                        if (dicOfCurve.ContainsKey(cPep))
                            elution.Add(cPep, dicOfCurve[cPep].eCurvePerMs.Area);
                        else
                            elution.Add(cPep, 0);
                    double tmp = 1.0;
                    if (elution.Count > 1)
                        tmp = Math.Abs(MathNet.Numerics.Statistics.Correlation.Pearson(average.Values, elution.Values));

                    dicOfCorrelations[dicOfCurve] = tmp;
                    corrs.Add(tmp);
                }

                //Remove worst curves                
                if (corrs.Count > 1)
                {
                    Dictionary<Dictionary<Peptide, MaxFlowElutionCurve>, double> dicOfCurves2 = new Dictionary<Dictionary<Peptide, MaxFlowElutionCurve>, double>();

                    double medianCorr = MathNet.Numerics.Statistics.Statistics.Median(corrs);
                    double maxCorr = medianCorr + 0.5 * MathNet.Numerics.Statistics.Statistics.Variance(corrs);

                    foreach (Dictionary<Peptide, MaxFlowElutionCurve> dic in dicOfCurves.Keys)
                        if (dicOfCorrelations[dic] <= maxCorr)
                            dicOfCurves2.Add(dic, dicOfCurves[dic]);

                    lastDicOfCurves = dicOfCurves2;
                }
            }//End of While nbRun not exhausted

            Dictionary<Peptide, ElutionCurveMerger> cumulDic = new Dictionary<Peptide, ElutionCurveMerger>();
            foreach (Dictionary<Peptide, MaxFlowElutionCurve> dicOfCurve in lastDicOfCurves.Keys)
            {
                foreach (Peptide cPep in dicOfCurve.Keys)
                {
                    if (!cumulDic.ContainsKey(cPep))
                        cumulDic.Add(cPep, new ElutionCurveMerger());

                    cumulDic[cPep].AddCurve(dicOfCurve[cPep], dicOfCorrelations[dicOfCurve]);
                }
            }
            Dictionary<Peptide, MaxFlowElutionCurve> peptideRatiosNoSpike = new Dictionary<Peptide, MaxFlowElutionCurve>();
            foreach (Peptide cPep in cumulDic.Keys)
                peptideRatiosNoSpike.Add(cPep, cumulDic[cPep].Merge());

            return peptideRatiosNoSpike;
        }

        public class SolvedResult
        {
            public double Ratio;
            public double NbFitTimes;
            public SolvedResult(double ratio, double nbTimes)
            {
                this.Ratio = ratio;
                this.NbFitTimes = nbTimes;
            }
        }

        public static Dictionary<PeptideSpectrumMatch, SolvedResult> SolveFromFragmentScoreTheoMZ(int nbProductsToKeep,
                                            double precision, IEnumerable<MsMsPeak> capacity, MassTolerance tolerance,
                                            double PrecursorIntensityInCTrap,
                                            Query query,
                                            out double underFlow, out double percentError, IConSol ConSole)
        {
            List<Dictionary<double, double>> unitSpectrum = new List<Dictionary<double, double>>();
            foreach (PeptideSpectrumMatch psm in query.psms)
            {
                double[] arrayFrag = psm.ComputeAACoverage();
                Dictionary<double, double> individualSpectrum = new Dictionary<double, double>();
                foreach (ProductMatch match in psm.AllProductMatches)
                {
                    if (!individualSpectrum.ContainsKey(match.theoMz))
                        individualSpectrum.Add(match.theoMz, 0);

                    double addedFactor = individualSpectrum[match.theoMz];
                    if (match.Fragment == null)
                        addedFactor += arrayFrag[match.fragmentPos];
                    else if (match.Fragment.IsReverse)
                        for (int i = match.fragmentPos - 1; i < arrayFrag.Length; i++)
                            addedFactor += arrayFrag[i];
                    else
                        for (int i = match.fragmentPos - 1; i >= 0; i--)
                            addedFactor += arrayFrag[i];
                    individualSpectrum[match.theoMz] += addedFactor;
                }
                double sum = 0.0;
                foreach (double val in individualSpectrum.Values)
                    sum += val;
                foreach (double key in new List<double>(individualSpectrum.Keys))
                    individualSpectrum[key] /= sum;

                unitSpectrum.Add(individualSpectrum);
            }

            Dictionary<double, double> mixedSpectrum = new Dictionary<double, double>();
            foreach (Dictionary<double, double> dic in unitSpectrum)
                foreach (double key in dic.Keys)
                    if (!mixedSpectrum.ContainsKey(key))
                        mixedSpectrum.Add(key, precision);
            double nbKeys = (double)mixedSpectrum.Count;
            foreach (double key in new List<double>(mixedSpectrum.Keys))
                mixedSpectrum[key] /= nbKeys;
            //string seq = query.psms[0].Peptide.BaseSequence;
            //for (int i = 0; i < seq.Length; i++)
            //    mixedSpectrum.Add(i, precision);


            List<double> solution = new List<double>();
            double tmpUnderflow = 0;
            Utilities.Methods.GradientDescent.SolveMaxFlowStyle(unitSpectrum, mixedSpectrum, out solution, out tmpUnderflow, ConSole, 1);

            double sumOfIntensities = 0;
            foreach (double val in mixedSpectrum.Values)
                sumOfIntensities += val;

            underFlow = tmpUnderflow;
            List<SolvedResult> result = GetResultList(solution, precision, underFlow, sumOfIntensities);

            Dictionary<PeptideSpectrumMatch, SolvedResult> resultPerSample = new Dictionary<PeptideSpectrumMatch, SolvedResult>();
            int k = 0;
            foreach (PeptideSpectrumMatch key in query.psms)
            {
                resultPerSample.Add(key, result[k]);
                k++;
            }

            percentError = (underFlow / sumOfIntensities);
            return resultPerSample;
        }

        public static Dictionary<PeptideSpectrumMatch, SolvedResult> SolveFromFragmentScore(int nbProductsToKeep, 
                                            double precision, IEnumerable<MsMsPeak> capacity, MassTolerance tolerance, 
                                            double PrecursorIntensityInCTrap,
                                            Query query,
                                            out double underFlow, out double percentError, IConSol ConSole)
        {
            List<Dictionary<double, double>> unitSpectrum = new List<Dictionary<double, double>>();
            foreach(PeptideSpectrumMatch psm in query.psms)
            {
                double[] arrayFrag = psm.ComputeAACoverage();
                Dictionary<double, double> individualSpectrum = new Dictionary<double, double>();
                for (int i = 0; i < arrayFrag.Length; i++)
                    individualSpectrum.Add(i, arrayFrag[i]);
                unitSpectrum.Add(individualSpectrum);
            }

            Dictionary<double, double> mixedSpectrum = new Dictionary<double, double>();
            string seq = query.psms[0].Peptide.BaseSequence;
            for (int i = 0; i < seq.Length; i++)
                mixedSpectrum.Add(i, precision);


            List<double> solution = new List<double>();
            double tmpUnderflow = 0;
            Utilities.Methods.GradientDescent.SolveMaxFlowStyle(unitSpectrum, mixedSpectrum, out solution, out tmpUnderflow, ConSole, 1);

            double sumOfIntensities = 0;
            foreach (double val in mixedSpectrum.Values)
                sumOfIntensities += val;

            underFlow = tmpUnderflow;
            List<SolvedResult> result = GetResultList(solution, precision, underFlow, sumOfIntensities);

            Dictionary<PeptideSpectrumMatch, SolvedResult> resultPerSample = new Dictionary<PeptideSpectrumMatch, SolvedResult>();
            int k = 0;
            foreach (PeptideSpectrumMatch key in query.psms)
            {
                resultPerSample.Add(key, result[k]);
                k++;
            }

            percentError = (underFlow / sumOfIntensities);
            return resultPerSample;
        }

        public static Dictionary<CharacterizedPrecursor, SolvedResult> SolveFromSpectrumBKP(IEnumerable<CharacterizedPrecursor> ratiosToFit, int nbProductsToKeep,
                                            long precision, IEnumerable<MsMsPeak> capacity, MassTolerance tolerance,
                                            int returnType,//0 for max flow, 1 for best flow, 2 for average
                                            double PrecursorIntensityInCTrap,
                                            ref double overFlow, ref double underFlow, ref double errorInPercent, IConSol ConSole)
        {
            List<List<double>> solutions = new List<List<double>>();
            List<long> average = new List<long>();
            List<MsMsPeak> expandedCapacity = new List<MsMsPeak>();
            double sumOfProducts = 0;
            foreach (MsMsPeak peak in capacity)
            {
                double intensityNormed = peak.Intensity / PrecursorIntensityInCTrap;
                expandedCapacity.Add(new MsMsPeak(peak.MZ, intensityNormed * precision, peak.Charge));
                sumOfProducts += peak.Intensity;
            }
            List<List<ProductMatch>> tmpRatiosToFit = new List<List<ProductMatch>>();
            //foreach (List<ProductMatch> list in ratiosToFit.Values)
            foreach (CharacterizedPrecursor prec in ratiosToFit)
            {
                List<ProductMatch> pms = new List<ProductMatch>();
                foreach (ProductMatch pm in prec.Fragments[nbProductsToKeep])
                {
                    ProductMatch newPm = new ProductMatch(pm);
                    newPm.obsIntensity = newPm.normalizedIntensity;// *PrecursorIntensityInCTrap;
                    pms.Add(newPm);
                }
                tmpRatiosToFit.Add(pms);
            }

            double error = ComputeMaxFlow(tmpRatiosToFit, expandedCapacity, tolerance, ref solutions, ref errorInPercent, ref average, ConSole);

            double sumOfIntensities = 0;
            foreach (MsMsPeak peak in expandedCapacity)
                sumOfIntensities += peak.Intensity;

            overFlow = 0;
            underFlow = error;

            List<SolvedResult> result = null;
            switch (returnType)
            {
                case 0:
                    result = GetResultList(solutions[0], precision, underFlow, sumOfIntensities);
                    break;
                case 1:
                    result = GetResultList(solutions[1], precision, underFlow, sumOfIntensities);
                    break;
                case 2:
                    List<double> tmpAverage = new List<double>();
                    foreach (double val in average)
                        tmpAverage.Add(val);
                    result = GetResultList(tmpAverage, precision, underFlow, sumOfIntensities);
                    break;
            }
            Dictionary<CharacterizedPrecursor, SolvedResult> resultPerSample = new Dictionary<CharacterizedPrecursor, SolvedResult>();
            int i = 0;
            foreach (CharacterizedPrecursor key in ratiosToFit)
            {
                resultPerSample.Add(key, result[i]);
                i++;
            }
            return resultPerSample;
        }

        private static List<SolvedResult> GetResultList(List<double> solution, double precision, double underFlow, double sumOfIntensities)
        {
            List<SolvedResult> rez = new List<SolvedResult>();
            //double sumVal = 0.0;
            //foreach(double val in solution)
            //    sumVal += val;
            
            //sumVal += (underFlow / sumOfIntensities) * precision;//Counter intuitive, but according to test samples, it is less precise

            foreach (double val in solution)
                rez.Add(new SolvedResult(val / precision, val));

            return rez;
        }

        private static double ComputeMaxFlow(List<List<ProductMatch>> spikedMatches,
                                    List<MsMsPeak> mixedSpectrum, MassTolerance tolerance,
                                ref List<List<double>> optimalSolutions,
                                ref double percentError,
                                ref List<long> average, IConSol ConSole)
        {
            //Create dictionnary of usefull peaks
            Dictionary<float, double> mixedFragDic = new Dictionary<float, double>();
            foreach (List<ProductMatch> fragmentRatio in spikedMatches)
            {
                foreach (ProductMatch match in fragmentRatio)
                {
                    if (!mixedFragDic.ContainsKey((float)match.theoMz))
                    {
                        float closest = -1;
                        foreach (float key in mixedFragDic.Keys)
                            if (Math.Abs(Utilities.Numerics.CalculateMassError(match.theoMz, key, tolerance.Units)) <= tolerance.Value)
                                closest = key;
                        if (closest > 0)
                        {
                            //ConSole.WriteLine("Potential problem with selected fragment masses ");
                            match.theoMz = closest;
                        }
                        else
                            mixedFragDic.Add((float)match.theoMz, 0);
                    }
                }
            }

            //Fill dictionnary with Intensities
            List<float> keys = new List<float>(mixedFragDic.Keys);
            foreach (MsMsPeak peak in mixedSpectrum)
                foreach (float mz in keys)
                {
                    if (Math.Abs(Utilities.Numerics.CalculateMassError(peak.MZ, mz, tolerance.Units)) <= tolerance.Value)
                        mixedFragDic[mz] += peak.Intensity;
                }

            List<long> localFlows = new List<long>();
            foreach (List<ProductMatch> fragmentRatio in spikedMatches)
                localFlows.Add(FindLocalMaximumFlow(fragmentRatio, mixedFragDic));

            Dictionary<float, double> virtualSpectrum = BuildVirtualSpectrum(spikedMatches, localFlows, mixedFragDic);
            double overError = MaxFlowHelper.ComputeOverflow(virtualSpectrum, mixedFragDic);
            double underError = MaxFlowHelper.ComputeUnderflow(virtualSpectrum, mixedFragDic);
            double[] bestIndexes = new double[spikedMatches.Count];

            int iterSize = 1;
            double bestOverallError = double.MaxValue;
            List<long> bestLocalFlows = new List<long>();
            Random rnd = new Random();
            while (overError >= 1 && iterSize < 10000)//anything less than 1 is an acceptable solution
            {
                for (int index = 0; index < bestIndexes.Length; index++)
                    bestIndexes[index] = -1;

                for (int i = 0; i < spikedMatches.Count; i++)
                {
                    if (localFlows[i] > 0)
                    {
                        localFlows[i] -= iterSize;

                        virtualSpectrum = BuildVirtualSpectrum(spikedMatches, localFlows, mixedFragDic);
                        double tmpErrorMinus = MaxFlowHelper.ComputeUnderflow(virtualSpectrum, mixedFragDic);
                        double tmpErrorPlus = MaxFlowHelper.ComputeOverflow(virtualSpectrum, mixedFragDic);

                        double tmpFlowRate = Math.Abs(overError - tmpErrorPlus);
                        double underDiff = Math.Abs(underError - tmpErrorMinus);
                        if (underDiff >= 1)
                            tmpFlowRate /= underDiff;
                        bestIndexes[i] = tmpFlowRate;

                        localFlows[i] += iterSize;
                    }
                }

                //Pick pseudo randomly best index
                double worstFlowRate = 0.0;
                for (int index = 0; index < bestIndexes.Length; index++)
                    if (bestIndexes[index] > worstFlowRate)
                    {
                        worstFlowRate = bestIndexes[index];
                    }

                if (worstFlowRate > 0)
                {
                    int nbMatching = 0;
                    for (int index = 0; index < bestIndexes.Length; index++)
                        if (bestIndexes[index] >= worstFlowRate)
                            nbMatching++;

                    int iterChoice = rnd.Next(0, nbMatching - 1);
                    int iterNb = 0;
                    for (int index = 0; index < bestIndexes.Length; index++)
                        if (bestIndexes[index] >= worstFlowRate)
                        {
                            if (iterChoice == iterNb)
                                localFlows[index] -= iterSize;
                            iterNb++;
                        }
                    iterSize = 1;
                }
                else
                    iterSize++;

                virtualSpectrum = BuildVirtualSpectrum(spikedMatches, localFlows, mixedFragDic);
                overError = MaxFlowHelper.ComputeOverflow(virtualSpectrum, mixedFragDic);
                underError = MaxFlowHelper.ComputeUnderflow(virtualSpectrum, mixedFragDic);
                if (overError + underError < bestOverallError)
                {
                    bestLocalFlows = new List<long>(localFlows);
                    bestOverallError = overError + underError;
                }
            }//End of while overflow > 1
            optimalSolutions.Clear();

            List<double> newList = new List<double>();
            foreach (long localFlow in localFlows)
                newList.Add(localFlow);
            optimalSolutions.Add(newList);

            newList = new List<double>();
            foreach (long localFlow in bestLocalFlows)
                newList.Add(localFlow);
            optimalSolutions.Add(newList);

            //Compute average
            if (bestOverallError < double.MaxValue)
            {
                average.Clear();
                for (int i = 0; i < optimalSolutions[0].Count; i++)
                {
                    double sum = 0.0;
                    foreach (List<double> solution in optimalSolutions)
                        sum += solution[i];
                    double avg = sum / (double)optimalSolutions.Count;
                    average.Add((long)avg);
                }
            }

            //Compute expected error in percentage
            double sumOfIntensities = 0;
            foreach (double val in mixedFragDic.Values)
                sumOfIntensities += val;
            percentError = underError / sumOfIntensities;

            virtualSpectrum = BuildVirtualSpectrum(spikedMatches, localFlows, mixedFragDic);

            return MaxFlowHelper.ComputeUnderflow(virtualSpectrum, mixedFragDic);
        }

        private static Dictionary<float, double> BuildVirtualSpectrum(List<List<ProductMatch>> fragRatios, List<long> ratios, Dictionary<float, double> fragments)
        {
            Dictionary<float, double> virtualFrag = new Dictionary<float, double>(fragments.Count);
            for (int i = 0; i < ratios.Count; i++)
            {
                foreach (ProductMatch match in fragRatios[i])
                {
                    if (!virtualFrag.ContainsKey((float)match.theoMz))
                        virtualFrag.Add((float)match.theoMz, match.obsIntensity * ratios[i]);
                    else
                        virtualFrag[(float)match.theoMz] += match.obsIntensity * ratios[i];
                }
            }
            return virtualFrag;
        }

        private static long FindLocalMaximumFlow(List<ProductMatch> spikedMatches, Dictionary<float, double> fragments)
        {
            Dictionary<float, double> virtualFrag = new Dictionary<float, double>(fragments.Count);

            bool letsGo = false;
            double maxCumul = 0;
            foreach (ProductMatch match in spikedMatches)
            {
                if (!virtualFrag.ContainsKey((float)match.theoMz))
                    virtualFrag.Add((float)match.theoMz, match.obsIntensity);
                else
                    virtualFrag[(float)match.theoMz] += match.obsIntensity;
                if (virtualFrag[(float)match.theoMz] > 0)
                    letsGo = true;
                maxCumul += match.obsIntensity;
            }
            maxCumul *= 10000000;
            long nbCumul = 0;
            if (letsGo)
            {
                nbCumul = 1;
                while (nbCumul < maxCumul)
                {
                    foreach (float key in fragments.Keys)
                    {
                        if (virtualFrag.ContainsKey(key))
                        {
                            if (virtualFrag[key] * nbCumul > fragments[key] + 1)
                                return nbCumul;
                        }
                    }
                    nbCumul++;
                }
            }
            return nbCumul;
        }
    }
}
