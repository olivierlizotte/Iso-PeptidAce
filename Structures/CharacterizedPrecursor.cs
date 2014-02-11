using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using PeptidAce.Utilities;
using PeptidAce.Iso.Methods;

namespace PeptidAce.Iso.Structures
{
    /// <summary>
    /// Wraps info related to the synthetic peptide. Info extracted from the individual runs
    /// </summary>
    public class CharacterizedPrecursor : PrecursorIon
    {
        public Peptide Peptide;
        public PeptideSpectrumMatches Psms;
        //Dictionnaries with one entry per number of product
        public Dictionary<int, double> PrecursorLossNormalizeFactor;
        public Dictionary<int, ElutionCurve> FragmentNormalizor;        
        public Dictionary<int, List<ProductMatch>> Fragments;
        public List<ProductMatch> AllFragments;

        //List of fragment masses that disapear across the elution peak
        public List<double> UnstableFragmentMasses;

        //Dictionnary with the peptide spectrum matches normalization factor
        private Dictionary<PeptideSpectrumMatch, double> DicOfPsmFactor;

        public Dictionary<int, Dictionary<double, double>> NormalizedFragments;

        public CharacterizedPrecursor(Sample sample, DBOptions dbOptions, Peptide peptide, IEnumerable<Query> queries, double mz)
            : base(sample, queries, mz, -1)
        {
            this.Peptide = peptide;
            Psms = new PeptideSpectrumMatches();

            foreach (Query query in queries)
                if (query.sample == sample)
                {
                    Psms.Add(new PeptideSpectrumMatch(query, peptide, dbOptions));
                    Charge = query.precursor.Charge;
                }
            Psms.Sort(PeptideSpectrumMatches.AscendingRetentionTime);
            UnstableFragmentMasses = new List<double>();
            AllFragments = Psms.GetCombinedSpectrum(dbOptions, this.eCurveIntensityPerMS, null, ref UnstableFragmentMasses);
            AllFragments.Sort(ProductMatch.DescendingWeightComparison);

            DicOfPsmFactor = new Dictionary<PeptideSpectrumMatch, double>();
            foreach (PeptideSpectrumMatch psm in Psms)
                DicOfPsmFactor.Add(psm, psm.ProbabilityScore());

            Fragments = new Dictionary<int, List<ProductMatch>>();
            PrecursorLossNormalizeFactor = new Dictionary<int, double>();
            FragmentNormalizor = new Dictionary<int, ElutionCurve>();
        }

        /// <summary>
        /// Retrieves common and unique fragments based on precomputed list of characterized peptides
        /// </summary>
        /// <param name="precursors"></param>
        /// <param name="nbFragmentPerPrec"></param>
        /// <param name="dbOptions"></param>
        /// <returns></returns>
        private static Dictionary<double, int> GetCommonFragmentMz(IEnumerable<CharacterizedPrecursor> precursors, int nbFragmentPerPrec, DBOptions dbOptions)
        {
            List<double> unstableMasses = new List<double>();
            foreach (CharacterizedPrecursor cPrec in precursors)
                unstableMasses.AddRange(cPrec.UnstableFragmentMasses);

            Dictionary<double, int> DicOfFragments = new Dictionary<double, int>();
            foreach (CharacterizedPrecursor cPrec in precursors)
            {
                int cumulFrag = 0;
                for (int i = 0; cumulFrag < nbFragmentPerPrec && i < cPrec.AllFragments.Count; i++)
                {
                    bool keep = true;
                    foreach (double mass in unstableMasses)
                        if (PeptidAce.Utilities.Numerics.MzDifference(cPrec.AllFragments[i].theoMz, mass, dbOptions.productMassTolerance.Units) < dbOptions.productMassTolerance.Value)
                            keep = false;
                    if (keep)
                    {
                        cumulFrag++;
                        if (!DicOfFragments.ContainsKey(cPrec.AllFragments[i].theoMz))
                            DicOfFragments.Add(cPrec.AllFragments[i].theoMz, 1);
                        else
                            DicOfFragments[cPrec.AllFragments[i].theoMz]++;
                    }
                }
            }
            return DicOfFragments;
        }

        /// <summary>
        /// Returns true is a precursor was characterized with sufficient quality (for a given number of product)
        /// </summary>
        /// <param name="nbProductsToKeep"></param>
        /// <returns></returns>
        public bool IsValid(int nbProductsToKeep)
        {
            return Fragments.ContainsKey(nbProductsToKeep) && NormalizedFragments.ContainsKey(nbProductsToKeep) && PrecursorLossNormalizeFactor.ContainsKey(nbProductsToKeep);
        }

        /// <summary>
        /// Buids a list with the ProductMatch objects of this characterized precursor
        /// related to the masses that are common from other isomers
        /// </summary>
        /// <param name="dicOfCommonFragments"></param>
        /// <param name="dbOptions"></param>
        /// <returns></returns>
        private List<ProductMatch> GetCombinedMatches(Dictionary<double, int> dicOfCommonFragments, DBOptions dbOptions)
        {
            List<ProductMatch> matches = new List<ProductMatch>(dicOfCommonFragments.Count);

            foreach (double mz in dicOfCommonFragments.Keys)
            {
                bool found = false;
                foreach (ProductMatch match in AllFragments)
                    if (match.theoMz == mz)
                    {
                        matches.Add(new ProductMatch(match));
                        found = true;
                    }

                if (!found)
                {
                    double sumPsmFactor = 0;
                    ProductMatch newMatch = new ProductMatch();
                    newMatch.theoMz = mz;
                    newMatch.weight = 0;
                    newMatch.obsIntensity = 0;
                    newMatch.normalizedIntensity = 0;
                    foreach (PeptideSpectrumMatch psm in Psms)
                    {
                        double timePoint = psm.Query.spectrum.RetentionTimeInMin * 1000.0 * 60.0;
                        foreach (MsMsPeak peak in psm.Query.spectrum.Peaks)
                        {
                            if (Math.Abs(Utilities.Numerics.CalculateMassError(peak.MZ, mz, dbOptions.productMassTolerance.Units)) <= dbOptions.productMassTolerance.Value)
                            {
                                newMatch.obsIntensity += peak.Intensity * DicOfPsmFactor[psm];
                                newMatch.normalizedIntensity += (peak.Intensity / 
                                                                this.eCurveIntensityPerMS.GetLocalArea(timePoint, timePoint + psm.Query.spectrum.InjectionTime)) * DicOfPsmFactor[psm];
                                                                //psm.Query.spectrum.PrecursorIntensityPerMilliSecond * psm.Query.spectrum.InjectionTime)) * DicOfPsmFactor[psm];
                                sumPsmFactor += DicOfPsmFactor[psm];
                                newMatch.weight++;
                            }
                        }
                    }
                    if (newMatch.weight > 0)
                    {
                        newMatch.obsIntensity /= sumPsmFactor;
                        newMatch.normalizedIntensity /= sumPsmFactor;
                    }
                    newMatch.weight *= newMatch.normalizedIntensity;
                    matches.Add(newMatch);
                }
            }

            double averageNormedIntensity = 0.0;
            foreach (ProductMatch match in matches)
                averageNormedIntensity += match.normalizedIntensity;

            if (matches.Count > 0)
                averageNormedIntensity /= (double)matches.Count;

            //Keep only most intense fragments (5% of average normalized intensity)
            foreach (ProductMatch pm in matches)
                if (pm.normalizedIntensity < averageNormedIntensity * 0.1)//0.05
                {
                    pm.normalizedIntensity = 0;
                    pm.obsIntensity = 0;
                }

            return matches;
        }

        /// <summary>
        /// Computes a normalization factor amongs synthetic peptides. Assumes synthetic sample runs were injected
        /// the same amount of peptides.
        /// </summary>
        /// <param name="allCorrespondingPrec"></param>
        /// <param name="average"></param>
        /// <param name="keep"></param>
        /// <returns></returns>
        private double GetNormalizePrecursorFactor(IEnumerable<CharacterizedPrecursor> allCorrespondingPrec, out double average, out bool keep)
        {
            keep = false;
            double PrecursorLossNormalizeFactor = 1.0;
            average = eCurve.Area;

            //average = area * Norm => Norm = average/area
            if (eCurve.Area > 0)
            {
                //Normalize matches based on precursor intensity differences
                double cumulArea = 0.0;
                int nbNonZero = 0;
                foreach (CharacterizedPrecursor precursor in allCorrespondingPrec)
                {
                    if (precursor.eCurve.Area > 0)
                    {
                        nbNonZero++;
                        cumulArea += precursor.eCurve.Area;
                    }
                }
                if (nbNonZero > 0)
                {
                    keep = true;
                    average = cumulArea / (double)nbNonZero;
                    //PrecursorLossNormalizeFactor = Math.Log(average, 2) / Math.Log(this.eCurve.Area, 2);
                    PrecursorLossNormalizeFactor = Math.Log(average, 10) / Math.Log(this.eCurve.Area, 10);
                    //PrecursorLossNormalizeFactor = average / this.eCurve.Area;
                    //PrecursorLossNormalizeFactor = 1.0;
                }
            }

            if (PrecursorLossNormalizeFactor > 4) PrecursorLossNormalizeFactor = 4;
            if (PrecursorLossNormalizeFactor < 0.25) PrecursorLossNormalizeFactor = 0.25;

            return PrecursorLossNormalizeFactor;
        }

        /// <summary>
        /// Computes a curve [precursor observed intensity ; correction so that precursor computed intensity from fragments matches observed intensity]
        /// </summary>
        /// <param name="dbOptions"></param>
        /// <param name="nbProductsToKeep"></param>
        /// <returns></returns>
        private ElutionCurve GetNormalizingCurve(DBOptions dbOptions, int nbProductsToKeep)
        {
            List<CharacterizedPrecursor> Isomers = new List<CharacterizedPrecursor>();
            Isomers.Add(this);

            ElutionCurve curve = new ElutionCurve();
            foreach (Query query in this.Queries)
            {
                double timeInMilliSeconds = query.spectrum.RetentionTimeInMin * 60.0 * 1000.0;
                
                double underFlow = 0;
                double percentError = 0;
                double intInTrap = query.spectrum.PrecursorIntensityPerMilliSecond * query.spectrum.InjectionTime;
                Dictionary<CharacterizedPrecursor, PositionnalIsomerSolver.SolvedResult> finalRatios = PositionnalIsomerSolver.SolveFromSpectrum(Isomers, nbProductsToKeep, query.spectrum.Peaks, dbOptions.productMassTolerance,
                                                        intInTrap, query.spectrum.PrecursorIntensity,
                                                        out underFlow, out percentError, dbOptions.ConSole);

                if (percentError < 0.5 && finalRatios[this].NbFitTimes > 0)
                {
                    double ratio = intInTrap / finalRatios[this].NbFitTimes;
                    if (ratio > 0.5 && ratio < 2)
                        curve.AddPoint(intInTrap, intInTrap / finalRatios[this].NbFitTimes);
                }
            }
            curve.Compute(CurveType.LINEAR, true);
            return curve;
        }

        /// <summary>
        /// Computes the actual curve from data acquired through the AddPoint method
        /// </summary>
        /// <param name="isomers"></param>
        /// <param name="minNbProducts"></param>
        /// <param name="maxNbProducts"></param>
        /// <param name="dbOptions"></param>
        public static void Update(IEnumerable<CharacterizedPrecursor> isomers, int minNbProducts, int maxNbProducts, DBOptions dbOptions)//, double precision)
        {
            foreach (CharacterizedPrecursor prec in isomers)
                prec.NormalizedFragments = new Dictionary<int,Dictionary<double,double>>();

            for (int nbProduct = minNbProducts; nbProduct <= maxNbProducts; nbProduct++)
            {
                Dictionary<double, int> dicOfFrags = GetCommonFragmentMz(isomers, nbProduct, dbOptions);
                foreach (CharacterizedPrecursor prec in isomers)
                    prec.Fragments.Add(nbProduct, prec.GetCombinedMatches(dicOfFrags, dbOptions));
                                
                foreach (CharacterizedPrecursor prec in isomers)
                {
                    if (prec.Fragments.ContainsKey(nbProduct))
                    {
                        Dictionary<double, double> dic = new Dictionary<double, double>();
                        foreach (double key in dicOfFrags.Keys)
                        {
                            dic.Add(key, 0.0);
                            foreach (ProductMatch match in prec.Fragments[nbProduct])
                                if (Math.Abs(Utilities.Numerics.CalculateMassError(match.theoMz, key, dbOptions.productMassTolerance.Units)) <= dbOptions.productMassTolerance.Value)
                                    dic[key] += match.normalizedIntensity;
                        }
                        prec.NormalizedFragments.Add(nbProduct, dic);
                    }
                }

                foreach (CharacterizedPrecursor prec in isomers)
                {
                    if (!prec.NormalizeFragments(isomers, nbProduct, dbOptions, false, false))//, precision))//If normalization fails, ignore this product
                        prec.Fragments.Remove(nbProduct);
                }
            }
        }

        /// <summary>
        /// Creates a dictionnary of characterized precursors based on ProPheus assignement results
        /// </summary>
        /// <param name="spikedSamples"></param>
        /// <param name="spikedResult"></param>
        /// <param name="dbOptions"></param>
        /// <param name="nbMinFragments"></param>
        /// <param name="nbMaxFragments"></param>
        /// <returns></returns>
        public static Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> GetSpikedPrecursors(Samples spikedSamples, Result spikedResult, DBOptions dbOptions, int nbMinFragments, int nbMaxFragments)//, long precision)
        {
            Dictionary<double, Dictionary<Query, int>> mzKeys = new Dictionary<double, Dictionary<Query, int>>();
            foreach (Query query in spikedResult.queries)
            {
                foreach (PeptideSpectrumMatch psm in query.psms)
                {
                    double mz = Utilities.Numerics.MZFromMass(psm.Peptide.MonoisotopicMass, query.spectrum.PrecursorCharge);
                    if (!mzKeys.ContainsKey(mz))
                    {
                        bool found = false;
                        foreach (double key in mzKeys.Keys)
                            if (Math.Abs(Utilities.Numerics.CalculateMassError(mz, key, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                            {
                                mz = key;
                                found = true;
                            }
                        if (!found)
                            mzKeys.Add(mz, new Dictionary<Query, int>());
                    }
                    if (!mzKeys[mz].ContainsKey(query))
                        mzKeys[mz].Add(query, 1);
                }
            }

            Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> spikes = new Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>>();
            foreach (Sample spikedSample in spikedSamples)
            {
                Dictionary<double, PrecursorIon> DicOfSpectrumMasses = PrecursorIon.GetPrecursors(spikedResult, spikedSample, dbOptions, mzKeys.Keys);
                foreach (double mzKey in DicOfSpectrumMasses.Keys)
                {
                    if (mzKeys.ContainsKey(mzKey))
                    {
                        //Pick the best PSM for each sample/precursor pair
                        Dictionary<Peptide, double> DicOfProbabilityScores = new Dictionary<Peptide, double>();

                        foreach (Query query in mzKeys[mzKey].Keys)
                            if (query.sample == spikedSample)
                            {
                                foreach (PeptideSpectrumMatch psm in query.psms)
                                    if (!DicOfProbabilityScores.ContainsKey(psm.Peptide))
                                        DicOfProbabilityScores.Add(psm.Peptide, psm.ProbabilityScore());
                                    else
                                        DicOfProbabilityScores[psm.Peptide] += psm.ProbabilityScore();
                            }

                        Peptide bestPeptide = null;
                        double bestScore = double.MinValue;
                        foreach (Peptide keyPep in DicOfProbabilityScores.Keys)
                            if (DicOfProbabilityScores[keyPep] > bestScore)
                            {
                                bestScore = DicOfProbabilityScores[keyPep];
                                bestPeptide = keyPep;
                            }
                        if (bestPeptide != null)
                        {
                            CharacterizedPrecursor cPrec = new CharacterizedPrecursor(spikedSample, dbOptions, bestPeptide, mzKeys[mzKey].Keys, mzKey);
                            //Don't keep precursors if they are not well characterized (unfragmented or missasigned)
                            if (cPrec.AllFragments.Count >= cPrec.Peptide.Length - 2)
                            {
                                if (!spikes.ContainsKey(mzKey))
                                    spikes.Add(mzKey, new Dictionary<Sample, CharacterizedPrecursor>());
                                if (!spikes[mzKey].ContainsKey(spikedSample))
                                    spikes[mzKey].Add(spikedSample, cPrec);
                                else
                                    Console.WriteLine("Twice??");
                            }
                        }
                    }
                }//End of foreach mzKey
            }//End of foreach spiked sample

            //Normalize intensities based on average area of each precursor
            List<double> tmpKeys = new List<double>(spikes.Keys);
            foreach (double mzKey in tmpKeys)
            {
                if(spikes[mzKey].Count > 0)
                    CharacterizedPrecursor.Update(spikes[mzKey].Values, nbMinFragments, nbMaxFragments, dbOptions);//, precision);
                else
                    spikes.Remove(mzKey);
            }//*/
            return spikes;
        }

        /// <summary>
        /// Normalizes fragments and precursors based on synthetic peptide runs
        /// </summary>
        /// <param name="allCorrespondingPrec"></param>
        /// <param name="nbProductsToKeep"></param>
        /// <param name="dbOptions"></param>
        /// <param name="normalizePrecursor"></param>
        /// <param name="normalizeFragments"></param>
        /// <returns></returns>
        private bool NormalizeFragments(IEnumerable<CharacterizedPrecursor> allCorrespondingPrec, int nbProductsToKeep, DBOptions dbOptions, bool normalizePrecursor, bool normalizeFragments)//, double precision)
        {
            bool keepNbProds = true;
            if (eCurve.Area > 0)
            {
                if (normalizeFragments)
                {
                    List<double> keys = new List<double>(NormalizedFragments[nbProductsToKeep].Keys);
                    double area = this.eCurve.Area;

                    FragmentNormalizor.Add(nbProductsToKeep, GetNormalizingCurve(dbOptions, nbProductsToKeep));
                }

                if (normalizePrecursor)
                {
                    double average = 0;
                    PrecursorLossNormalizeFactor.Add(nbProductsToKeep, GetNormalizePrecursorFactor(allCorrespondingPrec, out average, out keepNbProds));                    
                }
                else
                    PrecursorLossNormalizeFactor.Add(nbProductsToKeep, 1.0);
            }
            return keepNbProds;
        }
    }
}
