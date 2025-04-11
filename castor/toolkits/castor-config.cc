#include "sOutputManager.hh"

using namespace std;

int main(int argc, char *argv[]) {
  cout << sOutputManager::GetInstance()->GetPathToConfigDir() << endl;
  return 0;
}
