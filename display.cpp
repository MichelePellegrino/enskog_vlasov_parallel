#include "display.hpp"

IOHandler::IOHandler
(DSMC* dsmc, std::istream& is, std::ostream& os):
Motherbase(dsmc), input(is), output(os)
{

}

IOHandler::IOHandler
(DSMC* dsmc):
Motherbase(dsmc), input(std::cin), output(std::cout)
{

}
