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

#ifdef DEBUG
const std::string SPLIT_SITE_WRONG_SIZE_ERROR = "In conversion from linear localPRG string to graph, splitting the string by the next var site resulted in the wrong number of intervals. Perhaps ordering of numbers in GFA is irregular?!";
#else
const std::string SPLIT_SITE_WRONG_SIZE_ERROR = "In conversion from linear localPRG string to graph, splitting the string by the next var site resulted in the wrong number of intervals. Perhaps ordering of numbers in GFA is irregular?!";
#endif

#ifdef DEBUG
const std::string SPLIT_SITE_ERROR = "In conversion from linear localPRG string to graph, splitting the string by the next var site resulted in the first interval being non alphabetic. Perhaps ordering of numbers in GFA is irregular?!";
#else
const std::string SPLIT_SITE_ERROR = "In conversion from linear localPRG string to graph, splitting the string by the next var site resulted in the first interval being non alphabetic. Perhaps ordering of numbers in GFA is irregular?!";
#endif

#ifdef DEBUG
const std::string INTERVAL_MISSING_FROM_PATH_ERROR = "Can't look up interval as it isn't a node of local PRG graph";
#else
const std::string INTERVAL_MISSING_FROM_PATH_ERROR = "Can't look up interval as it isn't a node of local PRG graph";
#endif
#endif
