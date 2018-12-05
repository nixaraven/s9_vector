# s9_vector
Add-on for [simongog/sdsl-lite](https://github.com/simongog/sdsl-lite) that defines a new compressed bitvector called `s9_vector`

Implements a new type of compressed bit vector with support of `rank1` and `select1` operations. 
`s9_vector` is based on a previous implementation that was made with [libcds](https://github.com/fclaude/libcds) in mind.

Converts the input `bit_vector` into a compressed bitvector using *gap encoding* and then compressing the resulting integer vector using *Simple 9* encoding.
Using blocks of sizes determined by a template parameter.


## Usage

First install [simongog/sdsl-lite](https://github.com/simongog/sdsl-lite) as usual. Then include the new bitvector header that is in `include/sdsl/s9_vector.hpp` as it's shown in the next example:

```C++
#include <iostream>
#include <sdsl/int_vector.hpp>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include "s9_vector.hpp"

using namespace std;
using namespace sdsl;

int main(){
    bit_vector b = bit_vector(80*(1<<20),0);
    for (size_t i=0; i<b.size(); i+=100)
        b[i] = 1;
    cout << "Tamaños en Bytes:\n" << endl;
    cout << "Original\t" << size_in_bytes(b) << endl;
   
    rrr_vector<63> rrrb(b);
    cout << "RRR\t\t" << size_in_bytes(rrrb) << endl;
    
    sd_vector<> sdb(b);
    cout << "SD\t\t" << size_in_bytes(sdb) << endl;

    s9_vector<> s9b(b);
    cout << "S9\t\t" << size_in_bytes(s9b) << endl;
    
    s9_vector<128, int_vector<32>> s9b32(b);
    cout << "S9(32)\t\t" << size_in_bytes(s9b32) << endl;

    return 0;
}
```
