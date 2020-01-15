// g++ -I/usr/local/include/openbabel-2.0/ -static  recall.cpp -o recall -L/usr/local/lib -lopenbabel -lz -linchi
#include <iostream>
#include <string>
#include <set>
#include <dirent.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>

void collect_inchi(std::set<std::string> &inchi1, const  std::string &name1, int &count)
{
  count = 0;
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
      count++;
      mol.Clear();
      notatend = obconversion.Read(&mol);
    }
}

template <class InputIterator1, class InputIterator2>
size_t size_intersection (InputIterator1 first1, InputIterator1 last1,
			  InputIterator2 first2, InputIterator2 last2)
{
  size_t result = 0;
  while (first1!=last1 && first2!=last2)
  {
    if (*first1<*first2) ++first1;
    else if (*first2<*first1) ++first2;
    else {
      ++result; ++first1; ++first2;
    }
  }
  return result;
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
	      std::set<std::string> inchi1,inchi2;
	      int count1, count2;
	      collect_inchi(inchi1,name1, count1);
	      collect_inchi(inchi2,name2, count2);
	      total += inchi1.size();
	      computed += count2 - (count1 - inchi1.size());
	      identical +=  size_intersection(inchi1.begin(), inchi1.end(), inchi2.begin(), inchi2.end());
	    }
	}
      closedir (dir);
    } 
  else 
    {
      std::cerr << "Unable to open directory " << argv[1] << std::endl;
      return 1;
    }

  // std::cout << total <<" "<< identical << " " << double(identical) / total << " " << double(identical) / computed << std::endl;

  return(0);
}
