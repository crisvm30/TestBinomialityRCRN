BinomialRCRNGraph:=module()

description "Testing unconditional binomiality of the steady state ideal of a reversible chemical reaction network using graphs";

option
    package
#object
    ;

# Methods:
export
    binomialViaGraph,
    createGraph,
	isBinomial
;

#Varibles:
# local
#     ;

uses StringTools,GraphTheory;

(* binomialViaGraph:
 * Procedure that receives a reversible chemical reaction network as a parameter and returns a string specifying if the steady state ideal of the network is unconditionally binomial or not.
 * The parameter 'network' is a list of lists where each sub list represents a reaction. Each reaction contains 2 strings that represent the reactant and product complexes of the reaction.
 * For example, the network: 2A->B,C->D will have the next representation:[["2A","B"],["C","D"]]
 *)
binomialViaGraph := proc(network::list)
	
	print("The network is:");
	print(network);
	local G := createGraph(network);
	local graphVertices := Vertices(G);
	local graphEdges := Edges(G);
	print("Number of Vertices");
	print(nops(graphVertices));
	print("Number of Edges");
	print(nops(graphEdges));
	
	local markedSpecies := table();
	for local actualRV in graphVertices do
		#If the actual vertex is a reaction vertex do
		if Search("⇋",actualRV) > 0 then

			#get the species vertices connected to the actual reaction
			local speciesFound := false;
			local currentSpeciesV := "";
			local connSpecies := Neighbors(G,actualRV);

			for local s in connSpecies do
				if not(assigned(markedSpecies[s])) then # if species vertex is not marked
					speciesFound := true;
					currentSpeciesV := s;
					break;
				end if;
			end do;
			if speciesFound then
				markedSpecies[currentSpeciesV] := "marked";

				#Get all reaction vertices connected with the marked species vertex
				local sReactions := Neighbors(G,currentSpeciesV);
				
				for local s2 in connSpecies do
					if s2 <> currentSpeciesV then
						local multX := -GetEdgeWeight(G, {s2, actualRV}) / GetEdgeWeight(G, {currentSpeciesV,actualRV});
						DeleteEdge(G, {s2, actualRV});	

						for local rv in sReactions do
							if rv <> actualRV then
								if HasEdge(G, {s2, rv}) then
									local coeff := GetEdgeWeight(G, {currentSpeciesV, rv})*multX + GetEdgeWeight(G, {s2, rv});
									if coeff <> 0 then
										SetEdgeWeight(G, {s2, rv}, coeff);
									else
										DeleteEdge(G, {s2, rv});
									end if;
								else
									local coeff := GetEdgeWeight(G, {currentSpeciesV, rv})*multX;
									AddEdge(G, [{s2, rv}, coeff]);
								end if;
							end if;
						end do;
					end if;
				end do;
			end if;			
		end if;
	end do;	
	local R := isBinomial(G);
	if R then
		return "Result: The steady state ideal of the network is unconditionally binomial";
	else
		return "Result: The steady state ideal of the network is NOT unconditionally binomial";
	end if;
end proc:

createGraph := proc(network::list)
	
	local G := Graph('weighted');
	local numOfReactions := nops(network);
	local speciesSet := {};
	local i,j,k,y,m;

	for i from 1 to numOfReactions do
		local reactionVertex := "";
		local speciesInReaction := table();
		for j from 1 to nops(network[i]) do
			
			#string with the current complex(without unnecessary blank spaces)
			local auxString := SubstituteAll(network[i][j]," ","");
			reactionVertex := cat(reactionVertex, auxString);
			if j <>  nops(network[i]) then
				reactionVertex := cat(reactionVertex, "⇋");
			end if;
			#list that will contain the objects of the complex
			local currSpeciesList := Split(auxString,"+");

			for k from 1 to nops(currSpeciesList) do
				y := 1;
				for m from 1 to length(currSpeciesList[k]) while IsDigit(currSpeciesList[k][m]) do
					y := m+1;
				end do;

				#Get the species and add it to the speciesSet and graph
				local species := substring(currSpeciesList[k], y .. length(currSpeciesList[k]));
				if not(evalb(species in speciesSet)) then
					speciesSet := speciesSet union {species};
					G := AddVertex(G, species);					
				end if;

				#Get the coefficient of the current species in the complex
				local coef;
				coef := 1;
				if y > 1 then
					coef := parse(substring(currSpeciesList[k], 1 .. (y-1)));
				end if;

				#If we are in the product complex then the coefficient should be negative
				if j = 2 then
					coef := -coef;
				end if;

				if not(assigned(speciesInReaction[species])) then 
					speciesInReaction[species] := coef;
				else
					speciesInReaction[species] := speciesInReaction[species]+coef;
				end if;


			end do;
		end do;
		
		G := AddVertex(G, reactionVertex);

		#Get all species in the current reactionVertex
		local auxSpeciesList := [indices(speciesInReaction)];
		for local s in auxSpeciesList do
		#if the species is in the same amount in both complexes, its coeff will be 0 and thus should not create an edge between vertices
			if speciesInReaction[s[1]] <> 0 then
				#add edge from the current reactionVertex to the expected species vertices
				AddEdge(G, [{reactionVertex,s[1]}, speciesInReaction[s[1]] ]);	
			end if;
			
		end do;			
	end do;					
	return G;
end proc:

isBinomial := proc(G::function)
	local result := true;
	local components := ConnectedComponents(G);
	for local c in components do
		if nops(c) > 2 then
			result := false;
			break;	
		end if;	
	end do;
	return result;
end proc:

#====END MODULE====
end module;

