BinomialRCRN:=module()

description "Testing unconditional binomiality of a reversible chemical reaction network";

option
    package
#object
    ;

# Methods:
export
    binomialViaRREF,
    createBinomialCoeffMatrix,
    computeODECoeff
;

#Varibles:
# local
#     ;

uses StringTools,LinearAlgebra;



(* binomialViaRREF:
 * Procedure that receives a reversible chemical reaction network as a parameter and returns a string specifying if the steady state ideal of the  network is unconditionally binomial or not.
 * The parameter 'network' is a list of lists where each sub list represents a reaction. Each reaction contains 2 strings that represent the reactant and product complexes of the reaction.
 * For example, the network: 2A->B,C->D will have the next representation:[["2A","B"],["C","D"]]
 *)
binomialViaRREF := proc(network::list)

	local m;
	m := createBinomialCoeffMatrix(network);
	print("The network is:");
	print(network);
	print("The binomial coefficient matrix of the network is:");
	print(m);
	
	local H:= ReducedRowEchelonForm(m);
	print("The Reduced Row Echelon Form of the matrix is:");
	print(H);
	
	local numRows := RowDimension(H);
	local numCols := ColumnDimension(H);

	local i,j,isBinomial;
	isBinomial := true;
	for i from 1 to numRows while isBinomial do
		local numOfEntries;
		numOfEntries := 0;
		for j from 1 to numCols while isBinomial do
			if H[i][j] <> 0 then
				numOfEntries := numOfEntries + 1;
				if numOfEntries > 1 then
					# If there is more than 1 entry in a row then it is not binomial
					isBinomial := false;
				end if;
			end if;
		end do;
	end do;

	if isBinomial then
		return "Result: The steady state ideal of the network is unconditionally binomial";
	else
		return "Result: The steady state ideal of the network is NOT unconditionally binomial";
	end if;
end proc:


(*createBinomialCoeffMatrix:
 * Procedure that receives a reversible chemical reaction network and returns a matrix that represents the binomial coefficient matrix  of the network.
 *)
createBinomialCoeffMatrix := proc(network::list)
	local numOfReactions := nops(network);
	local speciesSet,ODEcoeffList;
	speciesSet,ODEcoeffList := computeODECoeff(network);
	local listOfRows := [];

	local i,k,j;
	for i from 1 to nops(speciesSet) do
		local currentRow := [];
	
		currentRow := Array([seq(0, i = 1 .. numOfReactions)]);
	
		#Updating the entries in the rows with the corresponding coefficient
		for j from 1 to nops(ODEcoeffList[speciesSet[i]]) do
			local currReactionNum := ODEcoeffList[speciesSet[i]][j][1];
			currentRow[currReactionNum] :=  ODEcoeffList[speciesSet[i]][j][2];

		end do;
		currentRow := convert(currentRow,list);
		listOfRows := [op(listOfRows),currentRow];
	end do;
	local m:= Matrix(listOfRows);
	return m;
end proc:

(*computeODECoeff:
 * Procedure that receives a chemical reaction network and returns 2 objects: a set containing the species of the network, and a table containing the ODE coefficients for each species
 *)
computeODECoeff := proc(network::list)

	local numOfReactions := nops(network);
	local speciesSet := {};
	local ODEcoeffList := table();

	local i,j,k,y,m;
	for i from 1 to numOfReactions do
		for j from 1 to nops(network[i]) do

			#string with the current complex(without unnecessary blank spaces)
			local auxString := SubstituteAll(network[i][j]," ","");

			#list that will contain the objects of the complex
			local currSpeciesList := Split(auxString,"+");

			for k from 1 to nops(currSpeciesList) do
				y := 1;
				for m from 1 to length(currSpeciesList[k]) while IsDigit(currSpeciesList[k][m]) do
					y := m+1;
				end do;

				#Get the specie and add it to the speciesSet
				local specie := substring(currSpeciesList[k], y .. length(currSpeciesList[k]));
				speciesSet := speciesSet union {specie};

				#Get the coefficient of the current specie in the complex
				local coef;
				coef := 1;
				if y > 1 then
					coef := parse(substring(currSpeciesList[k], 1 .. (y-1)));
				end if;

				#If we are in the product complex then the coefficient should be negative
				if j = 2 then
					coef := -coef;
				end if;

				#Adding the coefficient of the current specie in the ODEcoeffList
				if type(ODEcoeffList[specie],'list') then
					local notFound := true;
					local b;
					for b from 1 to nops(ODEcoeffList[specie]) while notFound do
						if ODEcoeffList[specie][b][1] = i then
							local auxList :=  ODEcoeffList[specie][b];
							 ODEcoeffList[specie] := subsop(b=NULL,ODEcoeffList[specie]);
							ODEcoeffList[specie] := [op(ODEcoeffList[specie]),[auxList[1],auxList[2]+coef] ];
							notFound := false
						end if;
					end do;
					if notFound = true then
						ODEcoeffList[specie] := [op(ODEcoeffList[specie]),[i,coef] ]
					end if;
				else
					ODEcoeffList[specie] := [ [i,coef] ]
				end if;
			end do;
		end do;
	end do;

return speciesSet,ODEcoeffList;
end proc:


#====END MODULE====
end module;
