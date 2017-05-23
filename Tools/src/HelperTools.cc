#include "../interface/HelperTools.h"

HelperTools::HelperTools()
{
  
}

HelperTools::~HelperTools()
{
  
}


bool HelperTools::fexists(const char *filename)
{
  std::ifstream ifile(filename);
  return ifile.good();
}

std::string HelperTools::ConvertDoubleToString(double Number)
{
  std::ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  convert << Number;
  return convert.str();
}

std::string HelperTools::DotReplace(double var)
{
  std::string str = ConvertDoubleToString(var);
  std::replace(str.begin(), str.end(), '.', 'p');
  return str;
}

std::string HelperTools::ConvertIntToString(int Number)
{
  std::ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  convert << Number;
  return convert.str();
}

std::string HelperTools::ConvertIntToString(int Number, int pad)
{
  std::ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  if ( pad > 1 ) { convert << std::setw(pad) << std::setfill('0');}
  convert << Number;
  return convert.str();
}

std::string HelperTools::MakeTimeStamp()
{
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  
  int year = now->tm_year - 100;  /// + 1900 to get current year
  int month = now->tm_mon + 1;
  int day = now->tm_mday;
  int hour = now->tm_hour;
  int min = now->tm_min;
  //int sec = now->tm_sec;
  
  std::string year_str = ConvertIntToString(year, 2);
  std::string month_str = ConvertIntToString(month, 2);
  std::string day_str = ConvertIntToString(day, 2);
  std::string hour_str = ConvertIntToString(hour, 2);
  std::string min_str = ConvertIntToString(min, 2);
  //std::string sec_str = ConvertIntToString(sec, 2);
  
  std::string date_str = year_str + month_str + day_str + "_" + hour_str + min_str;
  return date_str;
}
