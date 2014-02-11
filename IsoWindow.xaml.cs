using FirstFloor.ModernUI.Windows.Controls;
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
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Reflection;

namespace PeptidAce.ModernUI
{/*
    class ProxyDomain : MarshalByRefObject
    {
        public Assembly GetAssembly(string AssemblyPath)
        {
            try
            {
                string folder = AppDomain.CurrentDomain.BaseDirectory +
                    "PreRequisites" + System.IO.Path.DirectorySeparatorChar;

                // Create application domain setup information
                AppDomainSetup domaininfo = new AppDomainSetup();
                domaininfo.ApplicationBase = folder;

                Dictionary<string, AssemblyName> dic = new Dictionary<string, AssemblyName>();
                AppDomain domain = AppDomain.CreateDomain("ProteoWizard", AppDomain.CurrentDomain.Evidence, domaininfo);
                foreach (AssemblyName refAsmName in Assembly.ReflectionOnlyLoadFrom(AssemblyPath).GetReferencedAssemblies())
                {
                    dic.Add(refAsmName.FullName, refAsmName);
                 //   domain.Load(refAsmName);
                }
                int sizeOfDic = 0;//dic.Count;
                while(sizeOfDic < dic.Count)
                {
                    sizeOfDic = dic.Count;
                    foreach (string assName in dic.Keys)
                    {
                        foreach (AssemblyName refAsmName in Assembly.ReflectionOnlyLoadFrom(assName).GetReferencedAssemblies())
                            if (!dic.ContainsKey(refAsmName.FullName))
                                dic.Add(refAsmName.FullName, refAsmName);
                    }
                }
                //foreach (AssemblyName assName in dic.Values)
                //    domain.Load(assName);

                return domain.Load(AssemblyName.GetAssemblyName(AssemblyPath));
                //return Assembly.LoadFrom(AssemblyPath);
            }
            catch (Exception ex)
            {
                throw new InvalidOperationException("Error", ex);
            }
        }
    }//*/

    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class IsoWindow : ModernWindow
    {
        List<System.Reflection.Assembly> Assemblies = new List<System.Reflection.Assembly>();
        public IsoWindow()
        {
            InitializeComponent();
        }
    }
}
