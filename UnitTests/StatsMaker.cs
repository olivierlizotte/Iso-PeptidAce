using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using PeptidAce.Iso.Methods;
using PeptidAce.Utilities;
using PeptidAce.Utilities.Methods;
using PeptidAce.Iso.Structures;

namespace PeptidAce.Iso.UnitTests
{
    public class StatsMaker
    {
        public static void ProjectMerge( PeptidAce.Utilities.Interfaces.IConSol console)
        {
            string strProjectAll    = @"C:\_IRIC\Data\NB\ProjectFile_EverythingReplicates_Oct.csv";      
            string project          = @"C:\_IRIC\Data\NB\ProjectTest_AllAce_Spiked_19Oct.csv";
            string fastaFile        = @"C:\_IRIC\Data\NB\peptide.fasta";
            DBOptions options       = PositionnalIsomerSolver.CreateOptions(fastaFile, @"C:\_IRIC\Data\NB\Units\", 8, 0.05, console);
            Samples samplesMixed    = new Samples(strProjectAll, 0, options);
            Samples samplesSynth    = new Samples(project, 0, options);
            
            PositionnalIsomerSolver newSolver = new PositionnalIsomerSolver();
            newSolver.precTolPpm = 15;
            newSolver.prodTolDa = 0.05;
            newSolver.nbMinFragments = 5;
            newSolver.nbMaxFragments = 5;

            string[] synths = new string[samplesSynth.Count];
            for (int i = 0; i < synths.Length; i++)
                synths[i] = samplesSynth[i].sSDF;

            string[] mixed = new string[samplesMixed.Count];
            for (int i = 0; i < mixed.Length; i++)
                mixed[i] = samplesMixed[i].sSDF;

            newSolver.Solve(synths, mixed, fastaFile, Utilities.vsCSV.GetFolder(mixed[0]), options.ConSole);

            //Precompute Spiked peptide identifications
            Result SpikedResult = Ace.Start(options, fastaFile, samplesSynth, false, false);

            Result mixedResult = Ace.Start(options, fastaFile, samplesMixed, false, false);

            //Compute all usable spiked peptides
            Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> characterizedPeptides = CharacterizedPrecursor.GetSpikedPrecursors(samplesSynth, SpikedResult, mixedResult, samplesMixed, options, newSolver.nbMinFragments, newSolver.nbMaxFragments);

            Dictionary<Sample, List<MixedPrecursor>> mixedPrecursors = new Dictionary<Sample, List<MixedPrecursor>>();
            foreach (Sample mixedSample in samplesMixed)
                mixedPrecursors.Add(mixedSample, MixedPrecursor.GetMixedPrecursors(mixedSample, mixedResult, options, characterizedPeptides));

            Dictionary<Sample, List<Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>>> results = new Dictionary<Sample, List<Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>>>();
            //Get the list of precursors to characterize
            foreach (Sample mixedSample in samplesMixed)
            {
                foreach (double keyMz in characterizedPeptides.Keys)
                {
                    //List<Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>> listOfRatios = new List<Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>>();
                    foreach (MixedPrecursor mPrec in mixedPrecursors[mixedSample])
                    {
                        if (mPrec.MZ == keyMz)
                        {
                            // Compute Max Flow for this precursor
                            Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve> ratios = PositionnalIsomerSolver.GetRatios(characterizedPeptides, mPrec, options, newSolver.nbMinFragments, newSolver.nbMaxFragments);
                            
                            if (!results.ContainsKey(mixedSample))
                                results.Add(mixedSample, new List<Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>>());
                            results[mixedSample].Add(ratios);
                        }
                    }
                }
            }

            List<CharacterizedPrecursor> precursors = new List<CharacterizedPrecursor>();
            foreach (Dictionary<Sample, CharacterizedPrecursor> dic in characterizedPeptides.Values)
                foreach (CharacterizedPrecursor cP in dic.Values)
                    precursors.Add(cP);
            
            //Create average of each characterized peptide plus standard deviance
            vsCSVWriter writerArea = new vsCSVWriter(@"C:\_IRIC\Data\NB\Merge\stats_Area.csv");

            string lineC = "Count,";
            foreach(CharacterizedPrecursor cP in precursors)
                lineC += cP.Peptide.Sequence + ",";
            lineC += "Intensity per ms,";
            foreach (CharacterizedPrecursor cP in precursors)
                lineC += cP.Peptide.Sequence + ",";
            lineC += "Standard Deviation Count,";
            foreach (CharacterizedPrecursor cP in precursors)
                lineC += cP.Peptide.Sequence + ",";
            lineC += "Standard Deviation per ms,";
            foreach (CharacterizedPrecursor cP in precursors)
                lineC += cP.Peptide.Sequence + ",";
            writerArea.AddLine(lineC);

            foreach(int cond in samplesMixed.GetConditions())
            {
                Dictionary<CharacterizedPrecursor, Dictionary<int, MaxFlowElutionCurve>> deconvoluted = new Dictionary<CharacterizedPrecursor, Dictionary<int, MaxFlowElutionCurve>>();
                string sampleName = "";
                foreach (Sample mixedSample in results.Keys)
                {
                    if(mixedSample.PROJECT.CONDITION == cond)
                    {
                        sampleName = vsCSV.GetFileName_NoExtension(mixedSample.sSDF);
                        foreach(Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve> ratio in results[mixedSample])
                        {
                            foreach (CharacterizedPrecursor cP in ratio.Keys)
                            {
                                if (ratio[cP].eCurveCount.Area > 0)
                                {
                                    if (!deconvoluted.ContainsKey(cP))
                                        deconvoluted.Add(cP, new Dictionary<int, MaxFlowElutionCurve>());

                                    if (deconvoluted[cP].ContainsKey(mixedSample.PROJECT.REPLICATE))
                                    {
                                        if (deconvoluted[cP][mixedSample.PROJECT.REPLICATE].eCurveCount.Area < ratio[cP].eCurveCount.Area)
                                            deconvoluted[cP][mixedSample.PROJECT.REPLICATE] = ratio[cP];
                                    }
                                    else
                                        deconvoluted[cP].Add(mixedSample.PROJECT.REPLICATE, ratio[cP]);
                                    //deconvoluted[cP].Add(ratio[cP]);
                                }
                            }
                        }
                    }
                }

                
                Dictionary<int, double> totalIntensityCount = new Dictionary<int,double>();
                Dictionary<int, double> totalIntensityPerMs = new Dictionary<int,double>();

                foreach (CharacterizedPrecursor cP in precursors)
                {
                    if (deconvoluted.ContainsKey(cP))
                    {
                        foreach(int keyRep in deconvoluted[cP].Keys)
                        //foreach (MaxFlowElutionCurve curve in deconvoluted[cP])
                        {
                            if (!totalIntensityCount.ContainsKey(keyRep))
                            {
                                totalIntensityCount.Add(keyRep, 0.0);
                                totalIntensityPerMs.Add(keyRep, 0.0);
                            }
                            MaxFlowElutionCurve curve = deconvoluted[cP][keyRep];
                            totalIntensityCount[keyRep] += curve.eCurveCount.Area;
                            totalIntensityPerMs[keyRep] += curve.eCurvePerMs.Area;
                        }
                    }
                }
                string lineArea = sampleName + ",";
                string lineMS = ",";
                string stdDevCount = "";
                string stdDevMS = "";
                //1) Compute an average out of the replicates
                foreach(CharacterizedPrecursor cP in precursors)                
                {
                    if (deconvoluted.ContainsKey(cP))
                    {
                        double averageAreaMS = 0;
                        double averageAreaCount = 0;
                        foreach (MaxFlowElutionCurve curve in deconvoluted[cP].Values)
                        {
                            averageAreaCount += curve.eCurveCount.Area;
                            averageAreaMS += curve.eCurvePerMs.Area;
                        }
                        if (averageAreaCount > 0)
                        {
                            averageAreaCount = (averageAreaCount / ((double)deconvoluted[cP].Count));
                            averageAreaMS = (averageAreaMS / ((double)deconvoluted[cP].Count));
                            
                            double deNormAverageCount = 0.0;
                            double deNormAveragePerMs = 0.0;
                            List<double> repAreaCount = new List<double>();
                            List<double> repAreaMS = new List<double>();
                            foreach(int keyRep in deconvoluted[cP].Keys)
                            {
                                MaxFlowElutionCurve curve = deconvoluted[cP][keyRep];
                                double tmpCount = (curve.eCurveCount.Area / totalIntensityCount[keyRep]) * averageAreaCount;
                                deNormAverageCount += tmpCount;
                                repAreaCount.Add(tmpCount);
                                double tmpPerMs = (curve.eCurvePerMs.Area / totalIntensityPerMs[keyRep]) * averageAreaMS;
                                deNormAveragePerMs += tmpPerMs;
                                repAreaMS.Add(tmpPerMs);
                            }

                            lineArea += (deNormAverageCount / ((double) repAreaCount.Count)) + ",";
                            lineMS += (deNormAveragePerMs / ((double) repAreaMS.Count)) + ",";
                            if (repAreaCount.Count > 1)
                            {
                                stdDevCount += MathNet.Numerics.Statistics.ArrayStatistics.StandardDeviation(repAreaCount.ToArray()) + ",";
                                stdDevMS += MathNet.Numerics.Statistics.ArrayStatistics.StandardDeviation(repAreaMS.ToArray()) + ",";
                            }
                            else
                            {
                                stdDevCount += ",";
                                stdDevMS += ",";
                            }
                        }
                        else
                        {
                            lineArea += ",";
                            lineMS += ",";
                            stdDevCount += ",";
                            stdDevMS += ",";
                        }
                    }
                    else
                    {
                        lineArea += ",";
                        lineMS += ",";
                        stdDevCount += ",";
                        stdDevMS += ",";
                    }
                }
                writerArea.AddLine(lineArea + lineMS + "," + stdDevCount + "," + stdDevMS);

                //2) Add replicates results (to use for standard deviation)
            }
            writerArea.WriteToFile();
        }
    }
}
