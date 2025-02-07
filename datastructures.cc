// Datastructures.cc
//
// Student name:
// Student email:
// Student number:

#include "datastructures.hh"

#include <list>
#include <random>
#include <cmath>
#include <stack>

std::minstd_rand rand_engine; // Reasonably quick pseudo-random generator

template <typename Type>
Type random_in_range(Type start, Type end)
{
    auto range = end-start;
    ++range;

    auto num = std::uniform_int_distribution<unsigned long int>(0, range-1)(rand_engine);

    return static_cast<Type>(start+num);
}

// Modify the code below to implement the functionality of the class.
// Also remove comments from the parameter names when you implement
// an operation (Commenting out parameter name prevents compiler from
// warning about unused parameters on operations you haven't yet implemented.)

Datastructures::Datastructures() : alphOrder_(false), distOrder_(false), pubsChanged_(true)
{
    all_affiliations_ = {};
    affIds_ = {};
    nameIds_ = {};
    distIds_ = {};
    nameAffiliationMap_ = {};
    coordAffiliationMap_ = {};
    distCoordMMap_ = {};
    all_publications_ = {};
    pubIds_ = {};
    allConnections_ = {};
    vec_visited_ = {};
}

Datastructures::~Datastructures()
{
    clear_all();
}

unsigned int Datastructures::get_affiliation_count()
{
    return all_affiliations_.size();
}

void Datastructures::clear_all()
{
    for (const auto& aff : all_affiliations_)
        delete aff.second;

    for (const auto& pub : all_publications_)
        delete pub.second;

    all_affiliations_.clear();
    all_publications_.clear();
    affIds_.clear();
    pubIds_.clear();
    nameIds_.clear();
    distIds_.clear();
    allConnections_.clear();
    vec_visited_.clear();

    nameAffiliationMap_.clear();
    coordAffiliationMap_.clear();
    distCoordMMap_.clear();

    alphOrder_ = false;
    distOrder_ = false;

    pubsChanged_ = true;
    connectionsChanged_ = true;

    make_not_visited();
}

std::vector<AffiliationID> Datastructures::get_all_affiliations()
{
    return affIds_;
}

bool Datastructures::add_affiliation(AffiliationID id, const Name &name, Coord xy)
{
    if (!findId(id)) {
        auto affToAdd = new Affiliation{id, name, xy};
        affToAdd->distanceToOrigin_ = coordDistance(xy);

        affIds_.push_back(id);
        all_affiliations_[id] = affToAdd;
        nameAffiliationMap_[name] = id;
        coordAffiliationMap_[xy] = id;

        distCoordMMap_.insert({coordDistance(xy), xy});

        alphOrder_ = false;
        distOrder_ = false;

        return true; }

    return false;
}

Name Datastructures::get_affiliation_name(AffiliationID id)
{
    if (findId(id))
        return all_affiliations_[id]->name_;

    return NO_NAME;
}

Coord Datastructures::get_affiliation_coord(AffiliationID id)
{
    if (findId(id))
        return all_affiliations_[id]->coord_;

    return NO_COORD;
}

std::vector<AffiliationID> Datastructures::get_affiliations_alphabetically()
{
    if(!nameAffiliationMap_.empty() && alphOrder_ == false) {
        nameIds_.clear();
        nameIds_.reserve(nameAffiliationMap_.size());

        // nameAffiliationMap_ has names as keys in alphabetical order
        // so we can just push_back them into nameIds_
        for (const auto& aff : nameAffiliationMap_) {
            nameIds_.push_back(aff.second);
        }
        alphOrder_ = true;
    }
    return nameIds_;
}

std::vector<AffiliationID> Datastructures::get_affiliations_distance_increasing()
{
    if (!all_affiliations_.empty() && distOrder_ == false) {
        distIds_.clear();
        distIds_.reserve(all_affiliations_.size());

        // distCoordMMap_ has names as keys in alphabetical order
        // so we can just push_back them into distIds_ using coordAffiliationMap_
        for (const auto& aff : distCoordMMap_) {
            distIds_.push_back(coordAffiliationMap_[aff.second]);
        }
        distOrder_ = true;
    }
    return distIds_;
}

AffiliationID Datastructures::find_affiliation_with_coord(Coord xy)
{
    if (coordAffiliationMap_.find(xy) != coordAffiliationMap_.end())
        return coordAffiliationMap_[xy];

    else
        return NO_AFFILIATION;
}

bool Datastructures::change_affiliation_coord(AffiliationID id, Coord newcoord)
{
    if (findId(id)) {

        Coord oldCoord = all_affiliations_[id]->coord_;

        // Erase old coordinate and add new coordinate to coordAffiliationMap_
        coordAffiliationMap_.erase(oldCoord);
        coordAffiliationMap_[newcoord] = id;

        // Add new coordinate to all_affiliations_
        all_affiliations_[id]->coord_ = newcoord;
        all_affiliations_[id]->distanceToOrigin_ = coordDistance(newcoord);

        // Erase the old coordinate from distCoordMMap_
        auto range = distCoordMMap_.equal_range(coordDistance(oldCoord));
        for (auto i = range.first; i != range.second; ++i)
            if(i->second == oldCoord) {
                distCoordMMap_.erase(i);
                break;
            }
        // Insert new coordinate to distCoordMMap_
        distCoordMMap_.insert({coordDistance(newcoord), newcoord});

        // Erase the coordinate from distIds_ and add new coordinate to correct position
        distIds_.erase(std::remove(distIds_.begin(), distIds_.end(), id), distIds_.end());
        auto insertPosition = std::lower_bound(distIds_.begin(), distIds_.end(), id,
                                    [this](const AffiliationID& id1, const AffiliationID& id2) {
                        return compDist(all_affiliations_[id1],all_affiliations_[id2]); });

        distIds_.insert(insertPosition, id);
        connectionsChanged_ = true;

        return true;
    }
    return false;
}

bool Datastructures::add_publication(PublicationID id, const Name &name, Year year, const std::vector<AffiliationID> &affiliations)
{
    if (!findId(id)) {
        all_publications_.emplace(id, new Publication{id, name, year, affiliations});
        pubsChanged_ = true;

        for (const auto& aff : affiliations)
        {
            all_affiliations_[aff]->publications_.push_back(id);
        }

        update_connection(affiliations);
        return true;
    }
    return false;
}

std::vector<PublicationID> Datastructures::all_publications()
{
    if(pubsChanged_) {
        pubIds_.clear();
        pubIds_.reserve(all_publications_.size());
        for(const auto& pub : all_publications_) {
            pubIds_.push_back(pub.first);
        }
        pubsChanged_ = false;
    }
    return pubIds_;
}

Name Datastructures::get_publication_name(PublicationID id)
{
    if (findId(id))
        return all_publications_[id]->name_;

    return NO_NAME;
}

Year Datastructures::get_publication_year(PublicationID id)
{
    if (findId(id))
        return all_publications_[id]->year_;

    return NO_YEAR;
}

std::vector<AffiliationID> Datastructures::get_affiliations(PublicationID id)
{
    if (findId(id))
        return all_publications_[id]->affiliations_;

    return {NO_AFFILIATION};
}

bool Datastructures::add_reference(PublicationID id, PublicationID parentid)
{
    if (findId(id) && findId(parentid)) {
        all_publications_[id]->referenced_ = parentid;
        all_publications_[parentid]->referencing_.push_back(id);
        return true;
    }
    return false;
}

std::vector<PublicationID> Datastructures::get_direct_references(PublicationID id)
{
    if (findId(id))
        return all_publications_[id]->referencing_;

    return {NO_PUBLICATION};
}

bool Datastructures::add_affiliation_to_publication(AffiliationID affiliationid, PublicationID publicationid)
{
    if (findId(affiliationid) && findId(publicationid)) {
        all_affiliations_[affiliationid]->publications_.push_back(publicationid);
        all_publications_[publicationid]->affiliations_.push_back(affiliationid);

        update_connection(all_publications_[publicationid]->affiliations_);
        return true;
    }
    return false;
}

std::vector<PublicationID> Datastructures::get_publications(AffiliationID id)
{
    if (findId(id))
        return all_affiliations_[id]->publications_;

    return {NO_PUBLICATION};
}

PublicationID Datastructures::get_parent(PublicationID id)
{
    if (findId(id))
        return all_publications_[id]->referenced_;

    return NO_PUBLICATION;
}

std::vector<std::pair<Year, PublicationID> > Datastructures::get_publications_after(AffiliationID affiliationid, Year year)
{
    if (findId(affiliationid)) {
        auto pubs = all_affiliations_[affiliationid]->publications_;
        std::vector<std::pair<Year, PublicationID>> pubsAfter = {};
        for(const auto& id : pubs) {
            if (all_publications_[id]->year_ >= year) {
                pubsAfter.emplace_back(all_publications_[id]->year_, id);
            }
        }
        return pubsAfter;
    }
    return {{NO_YEAR, NO_PUBLICATION}};
}

std::vector<PublicationID> Datastructures::get_referenced_by_chain(PublicationID id)
{
    if (findId(id)) {
        std::vector<PublicationID> chain;
        recursiveRefByChain(id, chain);
        return chain;
    }
    return {NO_PUBLICATION};
}

std::vector<PublicationID> Datastructures::get_all_references(PublicationID id)
{
    if (findId(id)) {
        std::vector<PublicationID> refs;
        recursiveReferences(id, refs);
        return refs;
    }
    return {NO_PUBLICATION};
}

std::vector<AffiliationID> Datastructures::get_affiliations_closest_to(Coord xy)
{
    std::vector<AffiliationID> closestVec;
    closestVec.reserve(3);

    if (all_affiliations_.empty())
        return closestVec;

    // Calculate distances and sort affiliations based on distance from xy
    std::set<std::pair<float, AffiliationID>> distanceAffiliationPairs;
    for (const auto& aff : all_affiliations_) {
        float distance = coordDistance(aff.second->coord_, xy);

        if (distanceAffiliationPairs.size() != 3) // Add new pairs until 3 in total
            distanceAffiliationPairs.emplace(distance, aff.first);

        else {
            auto iter_last = distanceAffiliationPairs.end();
            iter_last--; // Iterator pointing to the last element of distanceAffiliationPairs
            if ((iter_last->first) > distance) {
                distanceAffiliationPairs.erase(iter_last); // Erase last element if smaller distance is found
                distanceAffiliationPairs.emplace(distance, aff.first); // Add new smaller distance to pair
            }
        }
    }
    // Add affiliation IDs to the vector to be returned
    for (const auto& [distance, affiliation] : distanceAffiliationPairs) {
           closestVec.emplace_back(affiliation);
    }
    return closestVec;
}

bool Datastructures::remove_affiliation(AffiliationID id)
{
    if (findId(id))
    {
        nameAffiliationMap_.erase(all_affiliations_[id]->name_);
        coordAffiliationMap_.erase(all_affiliations_[id]->coord_);
        distCoordMMap_.erase(all_affiliations_[id]->distanceToOrigin_);

        delete all_affiliations_[id];
        all_affiliations_.erase(id);

        if (!affIds_.empty())
            affIds_.erase(std::remove(affIds_.begin(), affIds_.end(), id), affIds_.end());

        if (!nameIds_.empty())
            nameIds_.erase(std::remove(nameIds_.begin(), nameIds_.end(), id), nameIds_.end());

        if (!distIds_.empty())
            distIds_.erase(std::remove(distIds_.begin(), distIds_.end(), id), distIds_.end());

        for(const auto& pub : all_publications_) {
            auto AffIter = find(pub.second->affiliations_.begin(), pub.second->affiliations_.end(), id);
            if(AffIter != pub.second->affiliations_.end())
            {
                pub.second->affiliations_.erase(AffIter);
                break; // Publication can only have one affiliation so loop can be ended here
            }
        }
        return true;
    }
    return false;
}

PublicationID Datastructures::get_closest_common_parent(PublicationID id1, PublicationID id2)
{
    if (get_parent(id1) == NO_PUBLICATION || get_parent(id2) == NO_PUBLICATION)
        return NO_PUBLICATION;

    std::vector<PublicationID> parents1_vec = get_referenced_by_chain(id1);
    std::vector<PublicationID> parents2_vec = get_referenced_by_chain(id2);

    // Create a set of parents2 for more efficient finding of parentID
    std::unordered_set<PublicationID> parents2_set(parents2_vec.begin(), parents2_vec.end());

    // Find the common parent while keeping the order from parents1_vec
    for (const auto& parent : parents1_vec) {
        if (parents2_set.count(parent) > 0)
            return parent;
    }
    return NO_PUBLICATION;
}

bool Datastructures::remove_publication(PublicationID publicationid)
{
    if (findId(publicationid)) {

        delete all_publications_[publicationid];
        all_publications_.erase(publicationid);
        pubIds_.erase(std::remove(pubIds_.begin(), pubIds_.end(), publicationid), pubIds_.end());

        for(const auto& pub : all_publications_) {
            // Erase PublicationID if it was a parent for another publication
            if(pub.second->referenced_ == publicationid)
                pub.second->referenced_ = NO_PUBLICATION;

            // Erase PublicationID if it was a reference for some other publication
            auto PubIter = find(pub.second->referencing_.begin(), pub.second->referencing_.end(), publicationid);
            if(PubIter != pub.second->referencing_.end()) {
                pub.second->referencing_.erase(PubIter);
            }
        }
        // Erase PublicationID from Affiliations' list of Publications
        for(const auto& aff : all_affiliations_) {
            auto AffIter = find(aff.second->publications_.begin(), aff.second->publications_.end(), publicationid);
            if(AffIter != aff.second->publications_.end()) {
                aff.second->publications_.erase(AffIter);

                // Make sure that connection weights are correct
                auto connections = all_affiliations_[aff.first]->connections_;
                for(auto& conn : connections) {
                    conn.weight = get_weight(conn.aff1, conn.aff2); // Calculate new connection weight
                    if (conn.weight == 0) {
                        // If connection weight goes to 0 it means that connection is gone and it must be removed.
                        connections.erase(std::remove(connections.begin(), connections.end(), conn), connections.end());
                    }
                }
                connectionsChanged_ = true;
            }
        }
        return true;
    }
    return false;
}

int Datastructures::coordDistance(Coord xy1, Coord xy2)
{
    int dx = xy1.x - xy2.x;
    int dy = xy1.y - xy2.y;
    return sqrt(dx * dx + dy * dy);
}

bool Datastructures::compDist(const Affiliation* aff1, const Affiliation* aff2)
{
    if (aff1->distanceToOrigin_ == aff2->distanceToOrigin_)
        return aff1->coord_.y < aff2->coord_.y;

    return aff1->distanceToOrigin_ < aff2->distanceToOrigin_;
}

template <typename IDType>
bool Datastructures::findId(IDType id)
{
    if constexpr (std::is_same<IDType, PublicationID>::value) {
        auto pub_iter = all_publications_.find(id);
        if (pub_iter != all_publications_.end()) {
            return true;
        }
    } else if constexpr (std::is_same<IDType, AffiliationID>::value) {
        auto aff_iter = all_affiliations_.find(id);
        if (aff_iter != all_affiliations_.end()) {
            return true;
        }
    }
    return false;
}

std::vector<PublicationID> Datastructures::recursiveRefByChain(PublicationID id, std::vector<PublicationID>& parents)
{
    if(all_publications_[id]->referenced_ == NO_PUBLICATION)
        return parents;

    auto parentId = get_parent(id);
    parents.push_back(parentId);
    return recursiveRefByChain(parentId, parents);
}

std::vector<PublicationID> Datastructures::recursiveReferences(PublicationID id, std::vector<PublicationID>& refs)
{
    if (!all_publications_[id]->referencing_.empty()) {
        auto referencing = get_direct_references(id);

        // Save publications that are referenced to a vector
        refs.insert(refs.end(),referencing.begin(), referencing.end());

        for (const auto& ref : referencing) {
            refs = recursiveReferences(ref, refs); // Recursively call each reference
        }
    }
    return refs;
}

Path Datastructures::get_connected_affiliations(AffiliationID id)
{
    if (findId(id))
        return all_affiliations_[id]->connections_;

    return {};
}

Path Datastructures::get_all_connections()
{
    if(connectionsChanged_) {
        allConnections_.clear();

        for(const auto& aff : all_affiliations_)
        {
            for(const auto& conn : aff.second->connections_)
            {
                if(aff.first == min(conn.aff1, conn.aff2))
                    allConnections_.push_back(conn);
            }
        }
        connectionsChanged_ = false;
    }
    return allConnections_;
}

Path Datastructures::get_any_path(AffiliationID source, AffiliationID target)
{
    if (findId(source) && findId(target)) {
        make_not_visited();
        return dfs(source, target, true);
    }
    return {};
}

Path Datastructures::get_path_with_least_affiliations(AffiliationID source, AffiliationID target)
{
    if (findId(source) && findId(target))
    {
        make_not_visited();
        std::list<Affiliation*> listConn;

        auto sourceAff = all_affiliations_[source];

        listConn.push_front(sourceAff);
        sourceAff->visited_ = "gray";
        vec_visited_.push_back(sourceAff);

        while (!listConn.empty())
        {
            Affiliation* current = listConn.front();
            listConn.pop_front();

            current->visited_ = "gray";
            vec_visited_.push_back(current);

            for(const auto &connection : all_affiliations_[current->id_]->connections_)
            {
                Affiliation* neighbor = all_affiliations_[connection.aff2];

                if (neighbor->visited_ == "white")
                {
                    neighbor->visited_ = "gray"; // Mark connection as visited

                    neighbor->path_ = current->path_;
                    neighbor->path_.push_back(connection);

                    vec_visited_.push_back(neighbor);
                    listConn.push_front(neighbor);

                    if(neighbor->id_ == target) {
                        return neighbor->path_;
                    }
                }
            }
            current->visited_ = "black";
        }
    }
    return {};
}

Path Datastructures::dfs(AffiliationID source, AffiliationID target, bool returnFirst)
{
    auto& startAffiliation = all_affiliations_[source];
    startAffiliation->visited_ = "gray";
    vec_visited_.push_back(startAffiliation);
    Path sfPath = {};

    std::stack<std::tuple<AffiliationID, Path>> stack;
    stack.push({source, {}});

    float bestFric = std::numeric_limits<float>::max();

    while (!stack.empty())
    {
        auto [currentID, currentPath] = stack.top();
        stack.pop();

        auto& currentAffiliation = all_affiliations_[currentID];

        for (const auto& connection : currentAffiliation->connections_)
        {
            auto& neighbor = all_affiliations_[connection.aff2];

            if (neighbor->visited_ == "white")
            {
                Path newPath = currentPath;
                newPath.push_back(connection);

                if (connection.aff2 == target)
                {
                    if (returnFirst)
                    {
                        return newPath;
                    }

                    float pathFric = get_max_friction(newPath);

                    if (pathFric < bestFric || sfPath.size() == 0)
                    {
                        if (pathFric == bestFric && newPath.size() > sfPath.size())
                            continue;

                        sfPath = newPath;
                        bestFric = pathFric;
                    }
                }
                else
                {
                    neighbor->visited_ = "gray";
                    vec_visited_.push_back(neighbor);
                    stack.push({connection.aff2, newPath});
                }
            }
        }
        currentAffiliation->visited_ = "black";
    }
    return sfPath;
}

Path Datastructures::get_path_of_least_friction(AffiliationID source, AffiliationID target)
{
    if (findId(source) && findId(target))
    {
        make_not_visited();
        Path path = dfs(source, target, false);
        return path;
    }
    return {};
}

float Datastructures::get_max_friction(Path path)
{
    if(path.empty())
        return std::numeric_limits<float>::max();

    float max = 0;
    for(const auto &conn : path)
        if (conn.friction() > max)
            max = conn.friction();
    return max;
}

PathWithDist Datastructures::get_shortest_path(AffiliationID source, AffiliationID target)
{
    if (findId(source) && findId(target))
    {
        make_not_visited();

        auto& startAffiliation = all_affiliations_[source];
        startAffiliation->visited_ = "gray";
        vec_visited_.push_back(startAffiliation);

        typedef std::pair<int, Affiliation*> pi;
        std::priority_queue<pi, std::vector<pi>, std::greater<pi>> pqPath;

        startAffiliation->distFromSource = 0;
        pqPath.push({0, startAffiliation});
        while (!pqPath.empty())
        {
            auto [currentDistance, currentAff] = pqPath.top();
            pqPath.pop();

            if(currentAff->id_ == target) {
                return createPath(startAffiliation, currentAff);
            }

            for(const auto& connection : currentAff->connections_)
            {
                auto& neighbor = all_affiliations_[connection.aff2];

                int newDistance = currentDistance + connection.distance;
                if(newDistance < neighbor->distFromSource) {
                    neighbor->pi = connection;
                    neighbor->distFromSource = newDistance;
                    pqPath.push({newDistance, neighbor});
                }

                if (neighbor->visited_ == "white") {
                    neighbor->visited_ = "gray";
                    vec_visited_.push_back(neighbor);

                }
            }
            currentAff->visited_ = "black";
        }
    }
    return {};
}

PathWithDist Datastructures::createPath(Affiliation* source, Affiliation* target)
{
    PathWithDist pathWithDist = {};
    Affiliation* prevAff = all_affiliations_[target->pi.aff2];
    while(prevAff != source)
    {
        pathWithDist.push_back({prevAff->pi, prevAff->pi.distance});
        prevAff = all_affiliations_[prevAff->pi.aff1];
    }
    std::reverse(pathWithDist.begin(), pathWithDist.end());
    return pathWithDist;
}

void Datastructures::update_connection(const std::vector<AffiliationID> affiliations)
{
    if(affiliations.size() >= 2)
    {
        auto index = affiliations.begin();

        for(;index != affiliations.end(); index++)
        {
            for(const auto& aff : affiliations) {

                if (aff != *index) // Make sure we don't add connection to itself
                {
                    // Check if the connection already exists
                    auto& connections = all_affiliations_[*index]->connections_;
                    auto existingConnection = std::find_if(connections.begin(), connections.end(),
                        [index, &aff](const Connection& existing) {
                            return (existing.aff1 == *index && existing.aff2 == aff) ||
                                   (existing.aff1 == aff && existing.aff2 == *index);
                        });

                    if (existingConnection != connections.end())
                    {
                        // Update the weight if the connection exists
                        existingConnection->weight = get_weight(*index, aff);
                    }
                    else
                    {
                        // Add the connection if it doesn't exist
                        auto distance = coordDistance(all_affiliations_[*index]->coord_, all_affiliations_[aff]->coord_);
                        Connection connection{ *index, aff, 1, distance};
                        connections.push_back(connection);
                    }
                }
            }
        }
        connectionsChanged_ = true;
    }
}

void Datastructures::make_not_visited()
{
    for(const auto& aff : vec_visited_)
    {
        aff->visited_ = "white";
        aff->path_.clear();
        aff->pi = NO_CONNECTION;
        aff->distFromSource = std::numeric_limits<int>::max();
    }
    vec_visited_.clear();
}

Weight Datastructures::get_weight(AffiliationID source, AffiliationID target)
{
    std::vector<PublicationID> pub_intersection;
    std::vector<PublicationID> source_pubs = all_affiliations_[source]->publications_;
    std::vector<PublicationID> target_pubs = all_affiliations_[target]->publications_;

    sort(source_pubs.begin(), source_pubs.end());
    sort(target_pubs.begin(), target_pubs.end());

    std::set_intersection(source_pubs.begin(), source_pubs.end(), target_pubs.begin(), target_pubs.end(),
                          std::back_inserter(pub_intersection));

    return pub_intersection.size();
}

