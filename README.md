C++ version of ntuple software

To run, have all .h and .cxx files in the same directory, then run:



 - `rootcling -v4 -f mydict.cxx  -rmf libmydict.rootmap -rml libmydict.so  LinkDef.h`



 - `g++ -shared -o libmydict.so mydict.cxx `\``root-config --cflags --libs`\`` -fPIC`



After the dictionary is created for vector vector, then compile



 - `g++ -g -o analyse_haa analyse_haa.cxx plotting.cxx objects.cxx mydict.cxx event.cxx filereader.cxx $(root-config --libs --cflags)`



And run



`./analyse_haa`
