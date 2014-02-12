using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using PeptidAce.Utilities;
using PeptidAce;

namespace PeptidAce.Iso.Structures
{
    /// <summary>
    /// A precursor Ion is a bunch of spectrum taken across an elution curve. It might be
    /// a single peptide or composed of coeluting isomers.
    /// </summary>
    public class PrecursorIon
    {
        //Elution curve based on intensity count
        public ElutionCurve eCurveIntensityCount;
        //Elution curve computed from intensities per millisecond
        public ElutionCurve eCurveIntensityPerMS;
        //Theoretical mz for this precursor
        public double MZ;
        //Charge
        public int Charge;
        //List of queries (one query per spectrum)
        public Queries Queries;
        //Sample file associated to this precursor
        public Sample Sample;

        public PrecursorIon(Sample sample, IEnumerable<Query> queries, double mz, int charge)
        {
            this.MZ = mz;
            this.Charge = charge;
            this.Queries = new Queries();
            foreach (Query query in queries)
                if (query.sample == sample)
                    this.Queries.Add(query);

            this.Queries.Sort(Query.AscendingRetentionTimeComparison);
            Dictionary<double, double> dicOfTimeInMsVsIntensityPerMs = new Dictionary<double, double>();
            Dictionary<double, double> dicOfTimeInMsVsIntensityCount = new Dictionary<double, double>();
            foreach (Query query in this.Queries)
            {
                double time = query.spectrum.RetentionTimeInMin * 60.0 * 1000.0;
                dicOfTimeInMsVsIntensityPerMs.Add(time, query.spectrum.PrecursorIntensityPerMilliSecond);
                dicOfTimeInMsVsIntensityCount.Add(time, query.spectrum.PrecursorIntensity);
            }
            this.eCurveIntensityCount = ElutionCurve.Create(dicOfTimeInMsVsIntensityCount);//dicOfTimeInMsVsIntensityPerMs);
            this.eCurveIntensityPerMS = ElutionCurve.Create(dicOfTimeInMsVsIntensityPerMs);
            this.Sample = sample;
        }
        
        /// <summary>
        /// Computes a list of precursor from a sample, based on Propheus identification results
        /// </summary>
        /// <param name="result"></param>
        /// <param name="sample"></param>
        /// <param name="dbOptions"></param>
        /// <param name="keys"></param>
        /// <returns></returns>
        public static Dictionary<double, PrecursorIon> GetPrecursors(Result result, Sample sample, DBOptions dbOptions, IEnumerable<double> keys)
        {
            Dictionary<double, PrecursorIon> DicOfSpectrumMasses = new Dictionary<double, PrecursorIon>();            
            foreach (Query query in result.queries)
            {
                if (query.sample == sample)
                {
                    double foundKey = query.spectrum.PrecursorMZ;

                    bool foundInKeys = false;
                    foreach (double key in keys)
                    {
                        double distance = Math.Abs(Utilities.Numerics.CalculateMassError(query.spectrum.PrecursorMZ, key, dbOptions.precursorMassTolerance.Units));
                        if (distance <= dbOptions.precursorMassTolerance.Value)
                        {
                            if (!foundInKeys || distance < Math.Abs(Utilities.Numerics.CalculateMassError(query.spectrum.PrecursorMZ, foundKey, dbOptions.precursorMassTolerance.Units)))
                                foundKey = key;
                            foundInKeys = true;
                        }
                    }
                    if (!foundInKeys)
                        foreach (double key in DicOfSpectrumMasses.Keys)
                            if (Math.Abs(Utilities.Numerics.CalculateMassError(query.spectrum.PrecursorMZ, key, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                                foundKey = key;

                    if (!DicOfSpectrumMasses.ContainsKey(foundKey))
                    {
                        PrecursorIon precIon = new PrecursorIon(sample, new Queries(dbOptions), foundKey, -1);
                        precIon.Queries.Add(query);
                        DicOfSpectrumMasses.Add(foundKey, precIon);
                    }
                    else
                        DicOfSpectrumMasses[foundKey].Queries.Add(query);
                }
            }

            //Split similar precursor mass not eluting at the same time

            //aussi://retirer le processus de clustering de pAce
            return DicOfSpectrumMasses;
        }

        //Spilt a precursor ion if the elution curve has a 'pause' (intensities drop to zero for a while)
        public IEnumerable<PrecursorIon> SplitBasedOnTime(DBOptions dbOptions)
        {
            if(Queries.Count > 0)
            {
                this.Queries.Sort(Query.AscendingRetentionTimeComparison);
                List<double> timePoints = new List<double>();
                for (int i = 1; i < Queries.Count; i++)
                    timePoints.Add(Queries[i].spectrum.RetentionTimeInMin - Queries[i-1].spectrum.RetentionTimeInMin);

                double variance = MathNet.Numerics.Statistics.Statistics.UpperQuartile(timePoints);
                Queries newQ = new Queries(dbOptions);
                newQ.Add(Queries[0]);
                for(int i = 1; i < Queries.Count; i++)
                {
                    if(timePoints[i-1] > 10 * variance)
                    {
                        yield return new PrecursorIon(Sample, newQ, MZ, Charge);
                        newQ.Clear();
                    }
                    newQ.Add(Queries[i]);
                }
                yield return new PrecursorIon(Sample, newQ, MZ, Charge);
            }
        }
    }
}
