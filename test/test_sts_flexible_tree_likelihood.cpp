#include "gtest/gtest.h"

#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/Nucleotide/HKY85.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>

#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>

#include "util.h"
#include "likelihood_vector.h"
#include "online_util.h"
#include "simple_flexible_tree_likelihood.h"
#include "beagle_flexible_tree_likelihood.h"

#include <libhmsbeagle/beagle.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>

namespace sts { namespace test { namespace beagle_flexible_tree_likelihood {

const bpp::DNA dna;
constexpr double TOLERANCE = 1e-4;

std::unique_ptr<bpp::TreeTemplate<bpp::Node>> treeOfPath(const std::string& newick_path)
{
    bpp::Newick newick_io;
    std::unique_ptr<bpp::Tree> tree(newick_io.read(newick_path));
    std::unique_ptr<bpp::TreeTemplate<bpp::Node>> tt(new bpp::TreeTemplate<bpp::Node>(*tree));
    // Avoid warnings from bio++ when branches are too short
    for(bpp::Node* node : tt->getNodes()){
    	if(node->hasDistanceToFather() && node->getDistanceToFather() < 1e-5)node->setDistanceToFather(1e-5);
    }
    return tt;
}

std::unique_ptr<bpp::SiteContainer> alignment_of_fasta_path(const std::string& fasta_path, const bpp::Alphabet& alphabet)
{
    std::ifstream aln_stream(fasta_path);
    return std::unique_ptr<bpp::SiteContainer>(sts::util::read_alignment(aln_stream, &alphabet));
}

void testKnownTree(std::string fasta_path,
                     std::string newick_path,
                     bpp::SubstitutionModel& model,
                     bpp::DiscreteDistribution& rate_dist)
{
    using bpp::Node;
    using bpp::Tree;
    using bpp::TreeTemplate;
    using namespace std;

    unique_ptr<bpp::TreeTemplate<Node>> tt = treeOfPath(newick_path);
    unique_ptr<bpp::SiteContainer> aln = alignment_of_fasta_path(fasta_path, dna);
    unique_ptr<bpp::SitePatterns> sp(new bpp::SitePatterns(aln.get()));

    std::vector<int> node_ids = tt->getNodesId();
    std::sort(node_ids.begin(), node_ids.end());
    for(size_t i = 1; i < node_ids.size(); i++) {
        ASSERT_EQ(node_ids[i], node_ids[i-1] + 1);
    }

	std::vector<string> names = aln->getSequencesNames();
	size_t nameCounter = names.size();
	for(bpp::Node* node : tt->getNodes()){
		if(node->isLeaf()){
			size_t pos = find(names.begin(), names.end(), node->getName()) - names.begin();
			node->setId(static_cast<int>(pos));
		}
		else {
			node->setId(static_cast<int>(nameCounter));
			nameCounter++;
		}
	}
    
    
    // BEAGLE
    sts::online::BeagleFlexibleTreeLikelihood beagle_calculator(*sp, model, rate_dist);
    beagle_calculator.initialize(model, rate_dist, *tt);
    
    // STS
    sts::online::SimpleFlexibleTreeLikelihood calculator(*sp, model, rate_dist);
    calculator.initialize(model, rate_dist, *tt);
    
    const double beagle_ll = beagle_calculator.calculateLogLikelihood();
    const double ll = calculator.calculateLogLikelihood();
//     const size_t llCalls = beagle_calculator.numberOfBeagleUpdateTransitionsCalls();

    // This should not cause any more beagle operations to be executed.
//     beagle_calculator.calculateLogLikelihood();
//     ASSERT_EQ(llCalls, beagle_calculator.numberOfBeagleUpdateTransitionsCalls());

    // Make dirty
//     beagle_calculator.invalidateAll();
//     const double beagle_ll_cached = beagle_calculator.calculateLogLikelihood();
//     const size_t llCalls2 = beagle_calculator.numberOfBeagleUpdateTransitionsCalls();
//     ASSERT_NEAR(beagle_ll, beagle_ll_cached, TOLERANCE);
//     ASSERT_EQ(llCalls * 2, llCalls2);

    // Bio++
    bpp::DRHomogeneousTreeLikelihood like(*tt, &model, &rate_dist, false, false);
    like.setData(*aln);
    like.initialize();
    const double bpp_ll = -like.getValue();

    ASSERT_NEAR(beagle_ll, bpp_ll, TOLERANCE);
    
    ASSERT_NEAR(beagle_ll, ll, TOLERANCE);

//     const int b = beagle_calculator.getDistalBuffer(tt->getRootNode());
//     ASSERT_NEAR(beagle_calculator.logLikelihood(b), beagle_ll, TOLERANCE);
}

TEST(STSFlexibleTreeLikelihood, ThirtyJCConstant)
{
    bpp::JCnuc model(&dna);
    bpp::ConstantRateDistribution rates;
    testKnownTree("data/thirty.ma", "data/thirty.tree", model, rates);
}

TEST(STSFlexibleTreeLikelihood, ThirtyJCGamma4)
{
    bpp::JCnuc model(&dna);
    bpp::GammaDiscreteRateDistribution rates(4, 0.234);
    testKnownTree("data/thirty.ma", "data/thirty.tree", model, rates);
}

TEST(STSFlexibleTreeLikelihood, ThirtyHKYGamma4)
{
    bpp::HKY85 model(&dna, 2.0, 0.4, 0.2, 0.15, 0.25);
    bpp::GammaDiscreteRateDistribution rates(4, 0.234);
    testKnownTree("data/thirty.ma", "data/thirty.tree", model, rates);
}

void testAttachmentLikelihood(const std::string& tree_path, const std::string& fasta_path,
                              const bpp::SubstitutionModel& model, const bpp::DiscreteDistribution& rates)
{
    using namespace bpp;
    using namespace sts::online;
    using std::string;
    using std::unique_ptr;
    using std::vector;

    unique_ptr<TreeTemplate<Node>> tree = treeOfPath(tree_path);
    unique_ptr<SiteContainer> aln = alignment_of_fasta_path(fasta_path, dna);
    unique_ptr<bpp::SitePatterns> sp(new bpp::SitePatterns(aln.get()));

    ASSERT_FALSE(tree->isMultifurcating());
    
	std::vector<string> names = aln->getSequencesNames();
	size_t nameCounter = names.size();
	for(bpp::Node* node : tree->getNodes()){
		if(node->isLeaf()){
			size_t pos = find(names.begin(), names.end(), node->getName()) - names.begin();
			node->setId(static_cast<int>(pos));
		}
		else {
			node->setId(static_cast<int>(nameCounter));
			nameCounter++;
		}
	}
	
	tree->getRootNode()->getSon(0)->setDistanceToFather(tree->getRootNode()->getSon(0)->getDistanceToFather()+tree->getRootNode()->getSon(1)->getDistanceToFather());
	tree->getRootNode()->getSon(1)->setDistanceToFather(0);
	
    sts::online::BeagleFlexibleTreeLikelihood fullCalculator(*sp, model, rates);
    fullCalculator.initialize(model, rates, *tree);
    const double fullLogLikelihood = fullCalculator.calculateLogLikelihood();

	//std::cout<<bpp::TreeTemplateTools::treeToParenthesis(*tree)<<std::endl;

    for(const string& leafName : tree->getLeavesNames()) {
        bpp::TreeTemplate<Node> tmpTree(*tree);
        Node* n = tmpTree.getNode(leafName);
//         tmpTree.newOutGroup(n);
		if(tmpTree.getRootNode()->getSon(1) == n || tmpTree.getRootNode()->getSon(0) == n) continue;
		if(tmpTree.getRootNode()->getSon(1)==n->getFather()) continue;

        const double pendant = n->getDistanceToFather();// + sibling(n)->getDistanceToFather();
        // Sibling becomes the root.
        // Insert on first child of the sibling.
        // When sibling becomes root, this edge has all of the length.
        const bpp::Node* insertEdge = sibling(n);
//         ASSERT_EQ(static_cast<size_t>(2), sib->getNumberOfSons()) << "Sibling must be bifurcating (dropped " << leafName << ")";
//         const bpp::Node* insertEdge = sib->getSon(0);

        const double origInsertLength = insertEdge->getDistanceToFather();
        const double origParentInsertLength = insertEdge->getFather()->getDistanceToFather();
        
//         const double distal = sib->getSon(1)->getDistanceToFather();
		const double distal = origInsertLength;
		const double proximal = origParentInsertLength;
        TreeTemplateTools::dropLeaf(tmpTree, leafName);
        
        sts::online::BeagleFlexibleTreeLikelihood calculator(*sp, model, rates);
        calculator.initialize(model, rates, tmpTree);

// size_t nameCounter = names.size();
// 	for(bpp::Node* node : tmpTree.getNodes()){
// 		if(node->isLeaf()){
// 			size_t pos = find(names.begin(), names.end(), node->getName()) - names.begin();
// 			node->setId(static_cast<int>(pos));
// 		}
// 		else {
// 			node->setId(static_cast<int>(nameCounter));
// 			nameCounter++;
// 		}
// 	}
	        //fullCalculator.calculateLogLikelihood();
        ASSERT_NEAR(insertEdge->getDistanceToFather(),
                    origInsertLength+origParentInsertLength,
                    TOLERANCE);

        const double attLike = calculator.calculateLogLikelihood(*insertEdge, leafName, pendant, distal, proximal);
        EXPECT_NEAR(fullLogLikelihood, attLike, TOLERANCE) << "removing " << leafName;
    }
}

void testAttachmentLikelihood2(const std::string& tree_path, const std::string& fasta_path,
                               bpp::SubstitutionModel& model,  bpp::DiscreteDistribution& rates)
{
    using namespace bpp;
    using namespace sts::online;
    using std::string;
    using std::unique_ptr;
    using std::vector;

    unique_ptr<TreeTemplate<Node>> tree = treeOfPath(tree_path);
    unique_ptr<SiteContainer> aln = alignment_of_fasta_path(fasta_path, dna);
    unique_ptr<bpp::SitePatterns> sp(new bpp::SitePatterns(aln.get()));

    
//     if(tree->isMultifurcating()){
//     	tree->newOutGroup(tree->getLeaves()[0]);  
//     	tree->resetNodesId();
//     }
    ASSERT_FALSE(tree->isMultifurcating());
    
	std::vector<string> names = aln->getSequencesNames();
	size_t nameCounter = names.size();
	for(bpp::Node* node : tree->getNodes()){
		if(node->isLeaf()){
			size_t pos = find(names.begin(), names.end(), node->getName()) - names.begin();
			node->setId(static_cast<int>(pos));
		}
		else {
			node->setId(static_cast<int>(nameCounter));
			nameCounter++;
		}
	}
	//std::cout<<tree->getRootNode()->getSon(0)->getName()<<std::endl;
	//std::cout<<"root "<<tree->getRootNode()->getSon(1)->getName()<<std::endl;
	
	tree->getRootNode()->getSon(0)->setDistanceToFather(tree->getRootNode()->getSon(0)->getDistanceToFather()+tree->getRootNode()->getSon(1)->getDistanceToFather());
	tree->getRootNode()->getSon(1)->setDistanceToFather(1e-6);
	
    sts::online::BeagleFlexibleTreeLikelihood fullCalculator(*sp, model, rates);
    fullCalculator.initialize(model, rates, *tree);
    const double fullLogLikelihood = fullCalculator.calculateLogLikelihood();

    for(const string& leafName : tree->getLeavesNames()) {
        bpp::TreeTemplate<Node> tmpTree(*tree);
        Node* n = tmpTree.getNode(leafName);
//         tmpTree.newOutGroup(n);
		if(tmpTree.getRootNode()->getSon(1) == n || tmpTree.getRootNode()->getSon(0) == n) continue;
		if(tmpTree.getRootNode()->getSon(1)==n->getFather()) continue;

        const double pendant = n->getDistanceToFather();

        TreeTemplateTools::dropLeaf(tmpTree, leafName);
        
        size_t nameCounter = names.size();
	for(bpp::Node* nn : tmpTree.getNodes()){
		if(nn->isLeaf()){
			size_t pos = find(names.begin(), names.end(), nn->getName()) - names.begin();
			nn->setId(static_cast<int>(pos));
		}
		else {
			nn->setId(static_cast<int>(nameCounter));
			nameCounter++;
		}
	}
        
        sts::online::BeagleFlexibleTreeLikelihood beagle_calculator(*sp, model, rates);
        beagle_calculator.initialize(model, rates, tmpTree);
        
        sts::online::SimpleFlexibleTreeLikelihood calculator(*sp, model, rates);
        calculator.initialize(model, rates, tmpTree);
        
        vector<Node*> nodes = tmpTree.getNodes();
        for( int i = 0; i < nodes.size(); i++){
	        Node* node = nodes[i];
        	if(!node->hasFather() || node->getFather()->getSon(1)==node) continue;
        	
        	bpp::TreeTemplate<Node> tmpTree2(tmpTree);
        	vector<Node*> nodes2 = tmpTree2.getNodes();
        	Node* theNode = nodes2[i];
        	
    		Node* new_node = new Node(8, "node"+std::to_string(tmpTree2.getNumberOfNodes()));
        	Node* new_leaf = new Node(1, leafName);
    		new_node->addSon(new_leaf);
    		new_leaf->setDistanceToFather(pendant);

		    Node* father = theNode->getFather();

		    // branch lengths
    		const double d = theNode->getDistanceToFather();

	    	size_t pos = father->getSonPosition(theNode);
    		father->setSon(pos, new_node);
    		new_node->addSon(theNode);

  		  // Attachment branch lengths
    		new_node->setDistanceToFather(d/2);
    		theNode->setDistanceToFather(d/2);
    		
    		std::vector<string> names = aln->getSequencesNames();
	
	nameCounter = names.size();
	for(bpp::Node* nn : tmpTree2.getNodes()){
		if(nn->isLeaf()){
			size_t pos = find(names.begin(), names.end(), nn->getName()) - names.begin();
			nn->setId(static_cast<int>(pos));
		}
		else {
			nn->setId(static_cast<int>(nameCounter));
			nameCounter++;
		}
	}
    		
        	sts::online::BeagleFlexibleTreeLikelihood beagle_calculator2(*sp, model, rates);
        	beagle_calculator2.initialize(model, rates, tmpTree2);
        	
        	sts::online::SimpleFlexibleTreeLikelihood calculator2(*sp, model, rates);
        	calculator2.initialize(model, rates, tmpTree2);
        	
        	const double beagle_lnl2 = beagle_calculator2.calculateLogLikelihood();
        	const double lnl2 = calculator2.calculateLogLikelihood();
        	
        	const double beagle_attLike = beagle_calculator.calculateLogLikelihood(*node, leafName, pendant, d/2, d/2);
        	const double attLike = calculator.calculateLogLikelihood(*node, leafName, pendant, d/2, d/2);
        	
        	EXPECT_NEAR(beagle_lnl2, beagle_attLike, TOLERANCE) << "removing " << leafName;
        	EXPECT_NEAR(lnl2, attLike, TOLERANCE) << "removing " << leafName;
        	EXPECT_NEAR(lnl2, beagle_lnl2, TOLERANCE) << "removing " << leafName;
        	
        	bpp::DRHomogeneousTreeLikelihood like(tmpTree2, &model, &rates, false, false);
		    like.setData(*aln);
	    	like.initialize();
	    	const double bpp_ll = -like.getValue();
		    EXPECT_NEAR(lnl2, bpp_ll, TOLERANCE) << "removing " << leafName;
        }
    }
}

void testAttachmentLikelihoodDeriv(const std::string& tree_path, const std::string& fasta_path,
                               bpp::SubstitutionModel& model, bpp::DiscreteDistribution& rates)
{
    using namespace bpp;
    using namespace sts::online;
    using std::string;
    using std::unique_ptr;
    using std::vector;

    unique_ptr<TreeTemplate<Node>> tree = treeOfPath(tree_path);
    unique_ptr<SiteContainer> aln = alignment_of_fasta_path(fasta_path, dna);
    unique_ptr<bpp::SitePatterns> sp(new bpp::SitePatterns(aln.get()));

    ASSERT_FALSE(tree->isMultifurcating());
    
	std::vector<string> names = aln->getSequencesNames();

    for(const string& leafName : tree->getLeavesNames()) {
        bpp::TreeTemplate<Node> tmpTree(*tree);
        Node* n = tmpTree.getNode(leafName);
//         tmpTree.newOutGroup(n);
		if(tmpTree.getRootNode()->getSon(1) == n || tmpTree.getRootNode()->getSon(0) == n) continue;
		if(tmpTree.getRootNode()->getSon(1)==n->getFather()) continue;

        const double pendant = n->getDistanceToFather();// + sibling(n)->getDistanceToFather();
        // Sibling becomes the root.
        // Insert on first child of the sibling.
        // When sibling becomes root, this edge has all of the length.
        const bpp::Node* insertEdge = sibling(n);

		bpp::DRHomogeneousTreeLikelihood like(*tree, &model, &rates, false, false);
		like.setData(*aln); 
		like.enableDerivatives(true);
	    like.initialize();
	    like.setParameterValue("BrLen" + std::to_string(insertEdge->getId()), insertEdge->getDistanceToFather());
	    const double dd1 = -(like.getFirstOrderDerivative("BrLen" + std::to_string(insertEdge->getId())));
	    const double dd2 = -(like.getSecondOrderDerivative("BrLen" + std::to_string(insertEdge->getId())));
	    
        const double origInsertLength = insertEdge->getDistanceToFather();
        const double origParentInsertLength = insertEdge->getFather()->getDistanceToFather();
        
//         const double distal = sib->getSon(1)->getDistanceToFather();
		const double distal = origInsertLength;
		const double proximal = origParentInsertLength;
        TreeTemplateTools::dropLeaf(tmpTree, leafName);
        
        sts::online::BeagleFlexibleTreeLikelihood beagle_calculator(*sp, model, rates);
        beagle_calculator.initialize(model, rates, tmpTree);
        
        sts::online::SimpleFlexibleTreeLikelihood calculator(*sp, model, rates);
        calculator.initialize(model, rates, tmpTree);

		size_t nameCounter = names.size();
		for(bpp::Node* node : tmpTree.getNodes()){
			if(node->isLeaf()){
				size_t pos = find(names.begin(), names.end(), node->getName()) - names.begin();
				node->setId(static_cast<int>(pos));
			}
			else {
				node->setId(static_cast<int>(nameCounter));
				nameCounter++;
			}
		}
		double d1;
		double d2;
		calculator.calculateLogLikelihood(*insertEdge, leafName, pendant, distal, proximal);
        calculator.calculateDistalDerivatives(*insertEdge, leafName, pendant, distal, proximal, &d1, &d2);
        
		double beagle_d1;
		double beagle_d2;
		beagle_calculator.calculateLogLikelihood(*insertEdge, leafName, pendant, distal, proximal);
        beagle_calculator.calculateDistalDerivatives(*insertEdge, leafName, pendant, distal, proximal, &beagle_d1, &beagle_d2);
        
        
	    	
        EXPECT_NEAR(d1, dd1, TOLERANCE) << "removing " << leafName;
        EXPECT_NEAR(d2, dd2, TOLERANCE) << "removing " << leafName;
        EXPECT_NEAR(d1, beagle_d1, TOLERANCE) << "removing " << leafName;
        EXPECT_NEAR(d2, beagle_d2, TOLERANCE) << "removing " << leafName;
        //std::cout <<d2<<" "<<dd2<<std::endl;
    }
}

// TEST(STSBeagleFlexibleTreeLikelihoodMidEdgeFive, JukesCantorConstant)
// {
//     bpp::JCnuc model(&dna);
//     bpp::ConstantRateDistribution rates;
//     testAttachmentLikelihood("data/5taxon/5taxon.tre", "data/5taxon/5taxon.fasta", model, rates);
// }

TEST(STSFlexibleTreeLikelihoodMidEdgeThirty, JukesCantorGamma6)
{
    bpp::JCnuc model(&dna);
    bpp::GammaDiscreteRateDistribution rates(6, 0.234);
    testAttachmentLikelihood("data/thirty.tree", "data/thirty.ma", model, rates);
}

TEST(STSFlexibleTreeLikelihoodDistalDerivThirty, JukesCantorConstant)
{
    bpp::JCnuc model(&dna);
    bpp::ConstantRateDistribution rates;
    testAttachmentLikelihoodDeriv("data/thirty.tree", "data/thirty.ma", model, rates);
}

TEST(STSFlexibleTreeLikelihoodAllMidEdgeThirty, JukesCantorConstant)
{
    bpp::JCnuc model(&dna);
    bpp::ConstantRateDistribution rates;
    testAttachmentLikelihood2("data/thirty.tree", "data/thirty.ma", model, rates);
}

TEST(STSFlexibleTreeLikelihoodMidEdgeThirty, JukesCantorGamma2)
{
    bpp::JCnuc model(&dna);
    bpp::GammaDiscreteRateDistribution rates(2, 0.234);
    testAttachmentLikelihood("data/thirty.tree", "data/thirty.ma", model, rates);
}

TEST(STSFlexibleTreeLikelihoodMidEdgeThirty, HKY85Constant)
{
    bpp::HKY85 model(&dna, 2.0, 0.25, 0.25, 0.3, 0.3);
    bpp::ConstantRateDistribution rates;
    testAttachmentLikelihood("data/thirty.tree", "data/thirty.ma", model, rates);
}

TEST(STSFlexibleTreeLikelihoodMidEdgeThirty, HKY85Gamma6)
{
    bpp::HKY85 model(&dna, 2.0, 0.4, 0.2, 0.15, 0.25);
    bpp::GammaDiscreteRateDistribution rates(6, 0.234);
    testAttachmentLikelihood("data/thirty.tree", "data/thirty.ma", model, rates);
}

TEST(STSFlexibleTreeLikelihoodMidEdge5taxon, HKY85Gamma6)
{
    bpp::HKY85 model(&dna, 2.0, 0.4, 0.2, 0.15, 0.25);
    bpp::GammaDiscreteRateDistribution rates(6, 0.234);
    testAttachmentLikelihood("data/5taxon/5taxon.tre", "data/5taxon/5taxon.fasta", model, rates);
}

TEST(STSFlexibleTreeLikelihoodMidEdge5taxon, JukesCantorConstant)
{
    bpp::JCnuc model(&dna);
    bpp::ConstantRateDistribution rates;
    testAttachmentLikelihood("data/5taxon/5taxon.tre", "data/5taxon/5taxon.fasta", model, rates);
}

}}} // namespaces
