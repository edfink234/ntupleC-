C++ version of ntuple software

To run, have all .h and .cxx files in the same directory, then run the following in the root command line:



gInterpreter->GenerateDictionary("vector<vector<string> >", "vector");



.L plotting.cxx



.L objects.cxx


.L event.cxx


.L filereader.cxx


.L analyse_haa.cxx


analyse_haa()






Alternatively, one can compile like a regular c++ program:




clang++ -g -o analyse_haa analyse_haa.cxx plotting.cxx objects.cxx event.cxx filereader.cxx $(root-config --libs --cflags)
