// g++ -I/usr/local/include/openbabel-2.0/ -static  detect.cpp -o recall -L/usr/local/lib -lopenbabel -lz -linchi
#include <iostream>
#include <string>
#include <set>
#include <dirent.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>

void collect_inchi(std::set<std::string> &inchi1, const  std::string &name1)
{
  OpenBabel::OBConversion obconversion;
  obconversion.SetInFormat("sdf");
  obconversion.SetOutFormat("inchi");
  obconversion.SetOptions("K", obconversion.OUTOPTIONS);
  OpenBabel::OBMol mol;
  bool notatend = obconversion.ReadFile(&mol, name1);
  while (notatend)
    {
      std::string inchi = obconversion.WriteString(&mol);
      if (!inchi.empty())
	inchi1.insert(inchi);
      mol.Clear();
      notatend = obconversion.Read(&mol);
    }
}

void print_errors(const std::set<std::string> &inchi1, const  std::string &name2)
{
  OpenBabel::OBConversion obconversion;
  obconversion.SetInFormat("sdf");
  obconversion.SetOutFormat("inchi");
  obconversion.SetOptions("K", obconversion.OUTOPTIONS);
  OpenBabel::OBMol mol;
  bool notatend = obconversion.ReadFile(&mol, name2);
  int i = 0;
  while (notatend)
    {
      std::string inchi = obconversion.WriteString(&mol);
      if (inchi.empty() || inchi1.find(inchi) == inchi1.end())
	{
	  std::cout << name2 << " " << i << std::endl;
	}
      mol.Clear();
      notatend = obconversion.Read(&mol);
      i++;
    }
}



int main(int argc,char **argv)
{

  if(argc<3)
    {
      std::cerr << "Usage: " << argv[0] <<" ground_truth/ computed/" << std::endl;
      return 1;
    }
  
  OpenBabel::obErrorLog.StopLogging();
 
  std::string folder1(argv[1]);
  std::string folder2(argv[2]);
  DIR *dir;
  struct dirent *ent;
  size_t total = 0, identical = 0, computed = 0;
  if ((dir = opendir (folder1.c_str())) != NULL) 
    {
      while ((ent = readdir (dir)) != NULL) 
	{
	  std::string name1 = folder1 + ent->d_name;
	  std::string name2 = folder2 + ent->d_name;
	  if (name1.size() > 4 && name1.substr(name1.size()-4) == ".sdf")
	    {
	      std::set<std::string> inchi1;
	      collect_inchi(inchi1,name1);
	      print_errors(inchi1,name2);
	    }
	}
      closedir (dir);
    } 
  else 
    {
      std::cerr << "Unable to open directory " << argv[1] << std::endl;
      return 1;
    }


  return(0);
}
