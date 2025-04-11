#include "JPetManager/JPetManager.h"
#include "ConvertEvents.h"

int main(int argc, const char** argv)
{
  try
  {
    JPetManager& manager = JPetManager::getManager();

    manager.registerTask<ConvertEvents>("ConvertEvents");
    manager.useTask("ConvertEvents", "cat.evt", "reco.unk.evt");

    manager.run(argc, argv);
  }
  catch (const std::exception& except)
  {
    std::cerr << "Unrecoverable error occured:" << except.what() << "Exiting the program!" << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
