using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using PeptidAce.Utilities;

namespace PeptidAce.Iso.Structures
{
    /// <summary>
    /// This class holds a mixed spectrum elution curve (spectrum and scan intensities).
    /// The class also stores results of the deconvolution process
    /// </summary>
    public class MixedPrecursor : PrecursorIon
    {
        //Deconvoluted isomers elution curves
        public Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve> PeptideRatios;
        //public Dictionary<Peptide, MaxFlowElutionCurve> PeptideRatiosNoSpike;
        public MixedPrecursor(Sample sample, IEnumerable<Query> queries, double mz)
            : base(sample, queries, mz, -1)
        {
        }

        public MixedPrecursor(Sample sample, PrecursorIon precursorIon, double mz)
            : base(sample, precursorIon.Queries, mz, -1)
        {
        }

        /// <summary>
        /// Computes a list of potential mixed precursors.
        /// Based on mass and elution time, spectrum are assembled and elution curves are computed for every precursor 
        /// that has a coverage satisfying the minimum necessary to undergo the deconvolution process.
        /// </summary>
        /// <param name="mixedSample"></param>
        /// <param name="mixedResult"></param>
        /// <param name="dbOptions"></param>
        /// <param name="charPeptides"></param>
        /// <returns></returns>
        public static List<MixedPrecursor> GetMixedPrecursors(Sample mixedSample, Result mixedResult, DBOptions dbOptions, Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> charPeptides)
        {
            Dictionary<double, PrecursorIon> DicOfSpectrumMasses = PrecursorIon.GetPrecursors(mixedResult, mixedSample, dbOptions, charPeptides.Keys);
            //Dictionary<double, MixedPrecursor> DicOfMixedPrecursor = new Dictionary<double, MixedPrecursor>();
            List<MixedPrecursor> listOfMixedPrec = new List<MixedPrecursor>();
            foreach (double key in DicOfSpectrumMasses.Keys)
            {
                if (charPeptides.ContainsKey(key))
                {
                    foreach (PrecursorIon precIon in DicOfSpectrumMasses[key].SplitBasedOnTime(dbOptions))
                    {
                        MixedPrecursor mixedPrecursor = new MixedPrecursor(mixedSample, precIon, key);
                        
                        //Don't try to characterize mixed precursors if there is less than five scans
                        if (mixedPrecursor.Queries.Count > 4)
                        {
                            foreach (Query q in mixedPrecursor.Queries)
                                q.spectrum.PrecursorIntensityPerMilliSecond = mixedPrecursor.eCurveIntensityPerMS.InterpolateIntensity(q.spectrum.RetentionTimeInMin * 1000.0 * 60.0);// + 0.5 * q.spectrum.InjectionTime);
                            listOfMixedPrec.Add(mixedPrecursor);
                        }
                            //DicOfMixedPrecursor.Add(key, mixedPrecursor);
                    }
                }
            }
            return listOfMixedPrec;
        }

        /// <summary>
        /// Based on deconvoluted curves from different set of products, this functions reports an average elution curve for each
        /// Each solution is weigthed based on the correlation with the (always evolving) average curve.
        /// </summary>
        /// <param name="dicOfCurveErrorsP"></param>
        /// <returns></returns>
        public Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve> ComputePeptideRatios(Dictionary<Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>, double> dicOfCurveErrorsP)
        {
            Dictionary<Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>, double> dicOfCorrelations = new Dictionary<Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>, double>();
            foreach (Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve> dicOfCurve in dicOfCurveErrorsP.Keys)
                if(dicOfCurve.Count > 0)
                    dicOfCorrelations.Add(dicOfCurve, 1.0 / (double)dicOfCurveErrorsP.Count);

            Dictionary<Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>, double> lastDicOfCurves = dicOfCurveErrorsP;
            /*
            //Purge worst curves
            Dictionary<Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>, double> dicOfCurves = new Dictionary<Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>, double>();
            if (dicOfCurveErrorsP.Count > 1)
            {
                double median = MathNet.Numerics.Statistics.Statistics.Median(dicOfCurveErrorsP.Values);
                double maxMed = median;// +0.5 * MathNet.Numerics.Statistics.Statistics.Variance(dicOfCurveErrorsP.Values);
                foreach (Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve> dic in dicOfCurveErrorsP.Keys)
                    if (dicOfCurveErrorsP[dic] <= maxMed)
                        dicOfCurves.Add(dic, dicOfCurveErrorsP[dic]);
            }
            else
                dicOfCurves = dicOfCurveErrorsP;
            
            lastDicOfCurves = dicOfCurves;
            int nbRun = 1;
            while (nbRun > 0)
            {
                nbRun--;
                //dicOfCurves = lastDicOfCurves;

                //Normalize already computed correlation factors for the remaning curves (sum must equal 1)
                double sumOfCorr = 0.0;
                foreach (Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve> dicOfCurve in dicOfCurves.Keys)
                    sumOfCorr += dicOfCorrelations[dicOfCurve];

                foreach (Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve> dicOfCurve in dicOfCurves.Keys)
                    dicOfCorrelations[dicOfCurve] /= sumOfCorr;

                //Compute average from weighted curves
                Dictionary<CharacterizedPrecursor, double> average = new Dictionary<CharacterizedPrecursor, double>();
                foreach (Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve> dicOfCurve in dicOfCurves.Keys)
                {
                    Dictionary<CharacterizedPrecursor, double> areas = GetAreas(dicOfCurve);
                    foreach (CharacterizedPrecursor cPep in areas.Keys)
                    {
                        if (!average.ContainsKey(cPep))
                            average.Add(cPep, 0);
                        average[cPep] += areas[cPep] * dicOfCorrelations[dicOfCurve];
                    }
                }

                //Compute correlation between average and curves
                List<double> corrs = new List<double>();
                foreach (Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve> dicOfCurve in dicOfCurves.Keys)
                {
                    Dictionary<CharacterizedPrecursor, double> elution = new Dictionary<CharacterizedPrecursor, double>();
                    foreach (CharacterizedPrecursor cPep in average.Keys)
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
                    Dictionary<Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>, double> dicOfCurves2 = new Dictionary<Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>, double>();

                    double medianCorr = MathNet.Numerics.Statistics.Statistics.Median(corrs);
                    double maxCorr = medianCorr + 0.5 * MathNet.Numerics.Statistics.Statistics.Variance(corrs);

                    foreach (Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve> dic in dicOfCurves.Keys)
                        if (dicOfCorrelations[dic] <= maxCorr)
                            dicOfCurves2.Add(dic, dicOfCurves[dic]);

                    lastDicOfCurves = dicOfCurves2;
                }
            }//End of While nbRun not exhausted
            //*/
            Dictionary<CharacterizedPrecursor, ElutionCurveMerger> cumulDic = new Dictionary<CharacterizedPrecursor, ElutionCurveMerger>();
            foreach (Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve> dicOfCurve in lastDicOfCurves.Keys)
            {
                foreach (CharacterizedPrecursor cPep in dicOfCurve.Keys)
                {
                    if (dicOfCorrelations.ContainsKey(dicOfCurve))
                    {
                        if (!cumulDic.ContainsKey(cPep))
                            cumulDic.Add(cPep, new ElutionCurveMerger());

                        cumulDic[cPep].AddCurve(dicOfCurve[cPep], dicOfCorrelations[dicOfCurve]);
                    }
                }
            }
            PeptideRatios = new Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>();
            foreach (CharacterizedPrecursor cPep in cumulDic.Keys)
                PeptideRatios.Add(cPep, cumulDic[cPep].Merge());

            return PeptideRatios;
        }

        /// <summary>
        /// Return area under the curve for each isomer (without using synthetic peptide runs)
        /// </summary>
        /// <param name="curves"></param>
        /// <returns></returns>
        public static Dictionary<Peptide, double> GetAreas(Dictionary<Peptide, MaxFlowElutionCurve> curves)
        {
            Dictionary<Peptide, double> newCurves = new Dictionary<Peptide, double>();
            try
            {
                foreach (Peptide cPep in curves.Keys)
                    newCurves.Add(cPep, curves[cPep].eCurvePerMs.Area);// * cPep.PrecursorLossNormalizeFactor[curves[cPep].nbProducts]);// * cPep.FragmentLossNormalizeFactor);
            }
            catch (Exception ex)
            {
                Console.WriteLine("Error in curve estimation :::");
                Console.WriteLine(ex.StackTrace);
            }
            return newCurves;
        }


        /// <summary>
        /// Retrieves area under the curve for each peptide isomer
        /// </summary>
        /// <param name="curves"></param>
        /// <returns></returns>
        public static Dictionary<CharacterizedPrecursor, double> GetAreas(Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve> curves)
        {
            Dictionary<CharacterizedPrecursor, double> newCurves = new Dictionary<CharacterizedPrecursor, double>();
            try
                {
                foreach (CharacterizedPrecursor cPep in curves.Keys)
                    newCurves.Add(cPep, curves[cPep].eCurvePerMs.Area * cPep.PrecursorLossNormalizeFactor[curves[cPep].nbProducts]);// * cPep.FragmentLossNormalizeFactor);
            }catch(Exception ex)
            {
                Console.WriteLine("Error in curve estimation");
                Console.WriteLine(ex.StackTrace);
            }
            return newCurves;
        }
    }
    
}
