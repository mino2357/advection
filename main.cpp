//
// advection
//
// u > 0
//

#include <iostream>
#include <vector>

int main(){
    auto u = std::vector<double>(100);
    
    for(unsigned int i = 0; i < u.size(); i++){
        std::cout << u[i] << std::endl;
    }
}