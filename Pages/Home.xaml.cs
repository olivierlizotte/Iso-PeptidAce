using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;
using FirstFloor.ModernUI.Presentation;
using PeptidAce.Iso.Structures;

namespace PeptidAce.ModernUI.Pages
{
    /// <summary>
    /// Interaction logic for Home.xaml
    /// </summary>
    public partial class Home : UserControl
    {
        private static Home _PTR;
        public Home()
        {
            InitializeComponent();
            _PTR = this;
        }
                        
        public static void UpdateLinks(IEnumerable<Sample> charSamples, Dictionary<Sample, List<MixedPrecursor>> mixedPrecursors)
        {//Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> characterizedPeptides
            while (_PTR.ModernTabList.Links.Count > 1)
                _PTR.ModernTabList.Links.RemoveAt(1);

            foreach(Sample sample in mixedPrecursors.Keys)
            //foreach(double mzKey in characterizedPeptides.Keys)
            {
              //  if (characterizedPeptides[mzKey].Count > 0)
                {
                    _PTR.ModernTabList.Links.Add(new Link
                    {
                        DisplayName = sample.Name,//((float)mzKey).ToString(),
                        Source = new Uri("/Content/Tabs.xaml?Sample=" + sample.Name, UriKind.Relative)
                    });
                }
            }

            //Spiked Peptides

            foreach (Sample sample in charSamples)
            //foreach(double mzKey in characterizedPeptides.Keys)
            {
                //  if (characterizedPeptides[mzKey].Count > 0)
                {
                    _PTR.ModernTabList.Links.Add(new Link
                    {
                        DisplayName = sample.Name,//((float)mzKey).ToString(),
                        Source = new Uri("/Content/TabsPeptides.xaml?Sample=" + sample.Name, UriKind.Relative)
                    });
                }
            }
            /*
            Dictionary<Sample, bool> doneSamples = new Dictionary<Sample,bool>();
            foreach (double mzKey in characterizedPeptides.Keys)
            {
                if (characterizedPeptides[mzKey].Count > 0)
                {
                    foreach (Sample sample in characterizedPeptides[mzKey].Keys)
                    {
                        if (!doneSamples.ContainsKey(sample))
                        {
                            _PTR.ModernTabList.Links.Add(new Link
                            {
                                DisplayName = "Synthetic Peptides Profiles",
                                Source = new Uri("/Content/ResultsPeptides.xaml?Sample=" + sample.Name, UriKind.Relative)
                            });
                            doneSamples.Add(sample, true);
                        }
                    }
                }
            }//*/
        }
    }
}
