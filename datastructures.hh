// Datastructures.hh
//
// Student name:
// Student email:
// Student number:

#ifndef DATASTRUCTURES_HH
#define DATASTRUCTURES_HH

#include <cmath>
#include <string>
#include <vector>
#include <tuple>
#include <utility>
#include <limits>
#include <functional>
#include <exception>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <queue>

// Types for IDs
using AffiliationID = std::string;
using PublicationID = unsigned long long int;
using Name = std::string;
using Year = unsigned short int;
using Weight = int;
struct Connection;
// Type for a distance (in arbitrary units)
using Distance = int;

using Path = std::vector<Connection>;
using PathWithDist = std::vector<std::pair<Connection,Distance>>;

// Return values for cases where required thing was not found
AffiliationID const NO_AFFILIATION = "---";
PublicationID const NO_PUBLICATION = -1;
Name const NO_NAME = "!NO_NAME!";
Year const NO_YEAR = -1;
Weight const NO_WEIGHT = -1;

// Return value for cases where integer values were not found
int const NO_VALUE = std::numeric_limits<int>::min();

// Type for a coordinate (x, y)
struct Coord
{
    int x = NO_VALUE;
    int y = NO_VALUE;
};


// Example: Defining == and hash function for Coord so that it can be used
// as key for std::unordered_map/set, if needed
inline bool operator==(Coord c1, Coord c2) { return c1.x == c2.x && c1.y == c2.y; }
inline bool operator!=(Coord c1, Coord c2) { return !(c1==c2); } // Not strictly necessary

struct CoordHash
{
    std::size_t operator()(Coord xy) const
    {
        auto hasher = std::hash<int>();
        auto xhash = hasher(xy.x);
        auto yhash = hasher(xy.y);
        // Combine hash values (magic!)
        return xhash ^ (yhash + 0x9e3779b9 + (xhash << 6) + (xhash >> 2));
    }
};

// Example: Defining < for Coord so that it can be used
// as key for std::map/set
inline bool operator<(Coord c1, Coord c2)
{
    if (c1.y < c2.y) { return true; }
    else if (c2.y < c1.y) { return false; }
    else { return c1.x < c2.x; }
}

// Return value for cases where coordinates were not found
Coord const NO_COORD = {NO_VALUE, NO_VALUE};

struct Connection
{
    AffiliationID aff1 = NO_AFFILIATION;
    AffiliationID aff2 = NO_AFFILIATION;
    Weight weight = NO_WEIGHT;
    Distance distance = NO_VALUE;

    float friction() const
    {
        if (weight != 0)
            return 1.0f / static_cast<float>(weight);
        else
            return std::numeric_limits<float>::infinity(); // Handle division by zero
    }

    bool operator==(const Connection& c1) const{
        return (aff1==c1.aff1) && (aff2==c1.aff2) && (weight==c1.weight);
    }

    bool operator!=(const Connection& c1) const
    {
        return !(*this == c1);
    }

    bool operator<(const Connection& c1) const
    {
        return weight < c1.weight;
    }
};
const Connection NO_CONNECTION{NO_AFFILIATION,NO_AFFILIATION,NO_WEIGHT};


// Return value for cases where Distance is unknown
Distance const NO_DISTANCE = NO_VALUE;

// This exception class is there just so that the user interface can notify
// about operations which are not (yet) implemented
class NotImplemented : public std::exception
{
public:
    NotImplemented() : msg_{} {}
    explicit NotImplemented(std::string const& msg) : msg_{msg + " not implemented"} {}

    virtual const char* what() const noexcept override
    {
        return msg_.c_str();
    }
private:
    std::string msg_;
};

// This is the class you are supposed to implement

class Datastructures
{
public:
    Datastructures();
    ~Datastructures();

    // Estimate of performance: O(1)
    // Short rationale for estimate: Constant because it only returns size of 'all_affiliations_'
    unsigned int get_affiliation_count();

    // Estimate of performance: O(n)
    // Short rationale for estimate: Both for-loops are linear and all "clear" functions are also linear
    void clear_all();

    // Estimate of performance: O(n)
    // Short rationale for estimate: Linear if affIds_ needs to be updated. Constant if affIds_ doesn't need to be updated.
    std::vector<AffiliationID> get_all_affiliations();

    // Estimate of performance: O(n)
    // Short rationale for estimate: "findId" function uses "find" function which is linear on worst case. Everything else is constant.
    bool add_affiliation(AffiliationID id, Name const& name, Coord xy);

    // Estimate of performance: O(n)
    // Short rationale for estimate: "findId" is linear on worst case. Accessing the name of the affiliation is constant.
    Name get_affiliation_name(AffiliationID id);

    // Estimate of performance: O(n)
    // Short rationale for estimate: "findId" is linear. Accessing the coordinations of the affiliation is constant.
    Coord get_affiliation_coord(AffiliationID id);


    // We recommend you implement the operations below only after implementing the ones above

    // Estimate of performance: O(n)
    // Short rationale for estimate: Linear because of for-loop.
    std::vector<AffiliationID> get_affiliations_alphabetically();

    // Estimate of performance: O(n)
    // Short rationale for estimate: Linear because of for-loop.
    std::vector<AffiliationID> get_affiliations_distance_increasing();

    // Estimate of performance: O(n)
    // Short rationale for estimate: std::find is linear.
    AffiliationID find_affiliation_with_coord(Coord xy);

    // Estimate of performance: O(n)
    // Short rationale for estimate: 'erase/remove' O(n), 'find' O(n), 'std::lower_bound' O(log(n)), 'insert' O(n), 'compDistOrigo' O(1),
    bool change_affiliation_coord(AffiliationID id, Coord newcoord);


    // We recommend you implement the operations below only after implementing the ones above

    // Estimate of performance: O(n)
    // Short rationale for estimate: "findId" is linear on worst case.
    bool add_publication(PublicationID id, Name const& name, Year year, const std::vector<AffiliationID> & affiliations);

    // Estimate of performance: O(n)
    // Short rationale for estimate: Linear because of the for-loop. Constant if pubIds_ is up to date.
    std::vector<PublicationID> all_publications();

    // Estimate of performance: O(n)
    // Short rationale for estimate: "findId" is linear on worst case.
    Name get_publication_name(PublicationID id);

    // Estimate of performance: O(n)
    // Short rationale for estimate: "findId" is linear on worst case.
    Year get_publication_year(PublicationID id);

    // Estimate of performance: O(n)
    // Short rationale for estimate: "findId" is linear on worst case.
    std::vector<AffiliationID> get_affiliations(PublicationID id);

    // Estimate of performance: O(n)
    // Short rationale for estimate: "findId" is linear on worst case.
    bool add_reference(PublicationID id, PublicationID parentid);

    // Estimate of performance: O(n)
    // Short rationale for estimate: "findId" is linear on worst case.
    std::vector<PublicationID> get_direct_references(PublicationID id);

    // Estimate of performance: O(n)
    // Short rationale for estimate: "findId" is linear on worst case.
    bool add_affiliation_to_publication(AffiliationID affiliationid, PublicationID publicationid);

    // Estimate of performance: O(n)
    // Short rationale for estimate: "findId" is linear on worst case.
    std::vector<PublicationID> get_publications(AffiliationID id);

    // Estimate of performance: O(n)
    // Short rationale for estimate: "findId" is linear on worst case.
    PublicationID get_parent(PublicationID id);

    // Estimate of performance: O(n)
    // Short rationale for estimate: linear because of the for-loop that iterates through all the publications of that Affiliation
    std::vector<std::pair<Year, PublicationID>> get_publications_after(AffiliationID affiliationid, Year year);

    // Estimate of performance: O(n)
    // Short rationale for estimate: findId" is linear on worst case. RecursiveReferences is logarithmic. So linear on worst case.
    std::vector<PublicationID> get_referenced_by_chain(PublicationID id);


    // Non-compulsory operations

    // Estimate of performance: O(n^2)
    // Short rationale for estimate: for-loop * recursive function that calls itself n times.
    std::vector<PublicationID> get_all_references(PublicationID id);

    // Estimate of performance: O(n^2)
    // Short rationale for estimate: linear for-loop * linear erase/emplace on worst case
    std::vector<AffiliationID> get_affiliations_closest_to(Coord xy);

    // Estimate of performance: O(n)
    // Short rationale for estimate: Searching all affiliations and publications has a worst-case scenario that is linear.
    bool remove_affiliation(AffiliationID id);

    // Estimate of performance: O(n^2)
    // Short rationale for estimate: Linear for-loop * linear unorderes_set.count() on worst case
    PublicationID get_closest_common_parent(PublicationID id1, PublicationID id2);

    // Estimate of performance: O(n^2)
    // Short rationale for estimate: Linear for-loop (n) * find algorithm for affiliations publications (m)
    bool remove_publication(PublicationID publicationid);

    // PRG 2 functions:

    // Estimate of performance: O(n)
    // Short rationale for estimate: "findId" is linear on worst case.
    std::vector<Connection> get_connected_affiliations(AffiliationID id);

    // Estimate of performance: O(n^2)
    // Short rationale for estimate: for-loop of (n) * for-loop (m)
    std::vector<Connection> get_all_connections();

    // Estimate of performance: O(n^2)
    // Short rationale for estimate: findId (n) + make_not_visited (n) + dfs (n^2)
    Path get_any_path(AffiliationID source, AffiliationID target);

    // PRG2 optional functions

    // Estimate of performance: O(n)
    // Short rationale for estimate: findId (n) + make_not_visited (n) + bfs (n)
    Path get_path_with_least_affiliations(AffiliationID source, AffiliationID target);

    // Estimate of performance: O(n * log(n)
    // Short rationale for estimate: findId (n) + make_not_visited (n) + dfs (n) + get_max_friction (n)
    Path get_path_of_least_friction(AffiliationID source, AffiliationID target);

    // Estimate of performance: O(n * log(n)
    // Short rationale for estimate: main-loop (n) * adding to priority queue (log n)
    PathWithDist get_shortest_path(AffiliationID source, AffiliationID target);


private:
    // Struct for a publication
    struct Publication
    {
        PublicationID id_ = NO_PUBLICATION;
        Name name_ = NO_NAME;
        Year year_ = NO_YEAR;
        std::vector<AffiliationID> affiliations_ = {};
        PublicationID referenced_ = NO_PUBLICATION;
        std::vector<PublicationID> referencing_ = {};
    };

    // Struct for an affiliation
    struct Affiliation
    {
        AffiliationID id_ = NO_AFFILIATION;
        Name name_ = NO_NAME;
        Coord coord_ = NO_COORD;
        std::vector<PublicationID> publications_ = {};
        int distanceToOrigin_ = NO_DISTANCE;

        std::vector<Connection> connections_ = {};
        std::string visited_ = "white";

        Path path_ = {};

        Connection pi = NO_CONNECTION;
        int distFromSource = std::numeric_limits<int>::max();
    };

    // Data Structures where information about affiliations and publications is stored
    using AffiliationMap = std::unordered_map<AffiliationID, Affiliation*>;
    using PublicationMap = std::unordered_map<PublicationID, Publication*>;
    AffiliationMap all_affiliations_;
    PublicationMap all_publications_;

    std::vector<AffiliationID> affIds_;        // Contains all Affiliation IDs
    std::vector<AffiliationID> nameIds_;
    std::vector<AffiliationID> distIds_;

    std::vector<PublicationID> pubIds_;         // Contains all Publication IDs

    Path allConnections_;
    std::vector<Affiliation*> vec_visited_;

    std::map<const Name, AffiliationID> nameAffiliationMap_;
    // Map that has affiliation coordinate as key, and a pointer to corresponding affiliation-struct as value
    std::unordered_map<Coord, AffiliationID, CoordHash> coordAffiliationMap_;

    bool alphOrder_; // Flag to see if affIds_ is not in alphabetical order by names
    bool distOrder_; // Flag to see if new Affiliations have been added and affDistances_ needs to be updated
    bool pubsChanged_; // Flag to see if new publications have been added and pubIds_ needs to be updated
    bool connectionsChanged_; // Flag to see if connections have been changed

    // Calculates distance to origo if only xy1 is given.
    // Or distance between two coordinates, if given xy1 and xy2.
    int coordDistance(Coord xy1, Coord xy2={0,0});

    // Comparison function for comparing which coordinate is closer to origo
    bool compDist(const Affiliation* aff1, const Affiliation* aff2);

    // Helper function that checks if AffiliationID or PublicationID can be found from data structure
    template <typename IDType>
    bool findId(IDType id);

    // Recursive functions needed for some commands
    std::vector<PublicationID> recursiveRefByChain(PublicationID id, std::vector<PublicationID>& parents);
    std::vector<PublicationID> recursiveReferences(PublicationID id, std::vector<PublicationID>& refs);

    void update_connection(const std::vector<AffiliationID> affiliations);
    void make_not_visited();
    Weight get_weight(AffiliationID source, AffiliationID target);
    float get_max_friction(Path path);

    Path dfs(AffiliationID source, AffiliationID target, bool returnFirst);
    PathWithDist createPath(Affiliation* source, Affiliation* target);

    using DistCoordMultimap = std::multimap<int, Coord>;
    DistCoordMultimap distCoordMMap_;

};

#endif // DATASTRUCTURES_HH
