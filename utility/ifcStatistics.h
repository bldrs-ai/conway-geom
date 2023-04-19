  #include <map>
  
  
  namespace webifc::statistics {
    extern bool exportObjs;// 			    = false;
    extern bool collectStats;//		    = false;
    extern bool verboseStats;//		    = false;
    extern bool shouldPrintTypeInfo;//	= false;

    extern std::map<uint32_t, uint32_t> uniqueTypeDefs;

    extern int printTypeInfo(const char *format, ...);
    extern void collectStatistics(uint32_t ifcLineType);
  }