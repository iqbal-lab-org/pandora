#ifndef __ERRORMESSAGES_H_INCLUDED__   // if errormessages.h hasn't been included yet...
#define __ERRORMESSAGES_H_INCLUDED__

#ifdef DEBUG
const std::string NODE_EXISTS_ERROR = "In conversion from linear localPRG string to graph, tried to add a node that already exists!";
#else
const std::string NODE_EXISTS_ERROR = "In conversion from linear localPRG string to graph, tried to add a node that already exists!";
#endif

#ifdef DEBUG
const std::string NODE_MISSING_ERROR = "In conversion from linear localPRG string to graph, tried to add edge between 2 nodes, one of which was missing!";
#else
const std::string NODE_MISSING_ERROR = "In conversion from linear localPRG string to graph, tried to add edge between 2 nodes, one of which was missing!";
#endif
#endif
