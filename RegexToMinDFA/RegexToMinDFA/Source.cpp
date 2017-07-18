#include<iostream>
#include<vector>
#include<stack>
#include<set>
#include<queue>
#include<map>
#include<string>
#include <chrono>

using namespace std;
using namespace std::chrono;

const int alphabetSize = 39;

struct NFAState
{
	vector<int> statesWithLetter[alphabetSize]; //states we go to with letters
	vector<int> statesWithEpsilon; //states we go to with epsilon transitions
	bool isFinal = false;
};

NFAState q; //empty non-final NFA state

struct DFAState
{
	int stateWithLetter[alphabetSize];
	bool isFinal;

	DFAState() {
		memset(stateWithLetter, -1, alphabetSize);
		isFinal = false;
	}
};

DFAState p; //empty non-final DFA state

void epsilonClosure(int state, set<int> &Eq, vector<NFAState> NFA)
{
	for (unsigned int i = 0; i < NFA[state].statesWithEpsilon.size(); ++i)
	{
		if (Eq.count(NFA[state].statesWithEpsilon[i]) == 0) //O(logn) TODO
		{
			Eq.insert(NFA[state].statesWithEpsilon[i]);
			epsilonClosure(NFA[state].statesWithEpsilon[i], Eq, NFA);
		}
	}
}

set<int> nextState(int l, set<int> Eq, vector<NFAState> NFA)
{
	set<int> result;
	for (std::set<int>::iterator it = Eq.begin(); it != Eq.end(); ++it)
	{
		for (unsigned int j = 0; j < NFA[*it].statesWithLetter[l].size(); ++j)
		{
			result.insert(NFA[*it].statesWithLetter[l][j]);
		}
	}

	return result;
}

vector<DFAState> determinize(vector<NFAState> NFA, int q0)
{
	vector<DFAState> DFA;
	int DFASize = 0; //number of DFA states
	int r = 0;
	bool finalIs;
	set<int> next;
	queue<set<int>> queueDFA;
	map<set<int>, int> mapDFA;
	set<int> Eq;

	mapDFA[Eq] = -1;
	Eq.insert(q0);
	epsilonClosure(q0, Eq, NFA);
	mapDFA[Eq] = r++;
	queueDFA.push(Eq);

	while (!queueDFA.empty())
	{
		DFA.push_back(p);
		Eq = queueDFA.front();
		finalIs = false;
		for (set<int>::iterator it = Eq.begin(); it != Eq.end(); ++it)
		{
			if (NFA[*it].isFinal == true)
			{
				finalIs = true;
			}
		}

		DFA[DFASize].isFinal = finalIs;
		for (int i = 0; i < alphabetSize; ++i)
		{
			Eq = queueDFA.front();
			next = nextState(i, Eq, NFA);
			Eq = next;
			for (set<int>::iterator it = Eq.begin(); it != Eq.end(); ++it)
			{
				epsilonClosure(*it, Eq, NFA);
			}

			if (mapDFA.count(Eq) == 0) //O(logn) TODO
			{
				DFA[DFASize].stateWithLetter[i] = r;
				mapDFA[Eq] = r++;
				queueDFA.push(Eq);
			}
			else
			{
				DFA[DFASize].stateWithLetter[i] = mapDFA.find(Eq)->second; //O(logn) TODO
			}

			next.clear(); //O(n) TODO
		}

		queueDFA.pop();
		DFASize++;
	}

	bool isTotal = true;

	for (int i = 0; i < DFASize; ++i)
	{
		for (int j = 0; j < alphabetSize; ++j)
		{
			if (DFA[i].stateWithLetter[j] == -1)
			{
				isTotal = false;
				break;
			}
		}
	}

	if (!isTotal)
	{
		for (int i = 0; i < DFASize; ++i)
		{
			for (int j = 0; j < alphabetSize; ++j)
			{
				if (DFA[i].stateWithLetter[j] == -1)
				{
					DFA[i].stateWithLetter[j] = DFASize;
				}
			}
		}

		DFA.push_back(p);
		for (int j = 0; j < alphabetSize; ++j)
		{
			DFA[DFASize].stateWithLetter[j] = DFASize;
		}
	}

	return DFA;
}

void print(vector<DFAState> minDFA)
{
	cout << endl;
	int br = 0;
	int printedLetters = 0;
	while (br < alphabetSize)
	{
		cout << "State\t|";
		char c;
		printedLetters = 0;
		for (int i = br; i < br + 3 && i < 26; ++i)
		{
			c = 'A' + i;
			cout << "\t" << c << "\t|";
			printedLetters++;
		}
		for (int i = br + printedLetters; i < br + 3 && i < 36; ++i)
		{
			cout << "\t" << i - 26 << "\t|";
			printedLetters++;
		}
		if (br + printedLetters < br + 3 && br + printedLetters == 36)
		{
			cout << "\t \t|";
			cout << "\t,\t|";
			cout << "\t$\t|";
		}
		cout << "\tFinal\t|" << endl;

		for (int i = 0; i<(int)minDFA.size(); i++)
		{
			cout << i << "\t|\t";
			for (int j = br; j < br + 3 && j<alphabetSize; ++j)
			{
				cout << minDFA[i].stateWithLetter[j] << "\t|\t";
			}

			if (minDFA[i].isFinal)
			{
				cout << "final\t|";
			}
			else
			{
				cout << "\t|";
			}
			cout << endl;
		}
		cout << endl;
		br += 3;
	}
}

void letter(int i, vector<NFAState> &NFA, int &NFASize, stack<int> &endpoints)
{
	NFA.push_back(q);
	NFA.push_back(q);
	NFA[NFASize].statesWithLetter[i].push_back(NFASize + 1);
	endpoints.push(NFASize);
	NFASize++;
	endpoints.push(NFASize);
	NFASize++;
}

void unite(vector<NFAState> &NFA, int &NFASize, stack<int> &endpoints)
{
	NFA.push_back(q);
	NFA.push_back(q);
	int q4 = endpoints.top(); endpoints.pop();
	int q3 = endpoints.top(); endpoints.pop();
	int q2 = endpoints.top(); endpoints.pop();
	int q1 = endpoints.top(); endpoints.pop();
	NFA[NFASize].statesWithEpsilon.push_back(q1);
	NFA[NFASize].statesWithEpsilon.push_back(q3);
	NFA[q2].statesWithEpsilon.push_back(NFASize + 1);
	NFA[q4].statesWithEpsilon.push_back(NFASize + 1);
	endpoints.push(NFASize);
	NFASize++;
	endpoints.push(NFASize);
	NFASize++;
}

void concatenate(vector<NFAState> &NFA, stack<int> &endpoints)
{
	int q4 = endpoints.top(); endpoints.pop();
	int q3 = endpoints.top(); endpoints.pop();
	int q2 = endpoints.top(); endpoints.pop();
	int q1 = endpoints.top(); endpoints.pop();
	NFA[q2].statesWithEpsilon.push_back(q3);
	endpoints.push(q1);
	endpoints.push(q4);
}

void kleeneStar(vector<NFAState> &NFA, int &NFASize, stack<int> &endpoints)
{
	NFA.push_back(q);
	NFA.push_back(q);
	int q2 = endpoints.top();
	endpoints.pop();
	int q1 = endpoints.top();
	endpoints.pop();
	NFA[NFASize].statesWithEpsilon.push_back(q1);
	NFA[NFASize].statesWithEpsilon.push_back(NFASize + 1);
	NFA[q2].statesWithEpsilon.push_back(q1);
	NFA[q2].statesWithEpsilon.push_back(NFASize + 1);
	endpoints.push(NFASize);
	NFASize++;
	endpoints.push(NFASize);
	NFASize++;
}

void copyNFA(vector<NFAState> &NFA, vector<NFAState> &copy, int state, set<int> &states, int &copySize)
{
	copy[state] = NFA[state];
	states.insert(state);
	copySize++;
	for (unsigned int i = 0; i < NFA[state].statesWithEpsilon.size(); ++i)
	{
		if (states.count(NFA[state].statesWithEpsilon[i]) == 0)
		{
			copyNFA(NFA, copy, NFA[state].statesWithEpsilon[i], states, copySize);
		}
	}

	for (int i = 0; i < alphabetSize; ++i)
	{
		for (int j = 0; j < NFA[state].statesWithLetter[i].size(); ++j)
		{
			if (states.count(NFA[state].statesWithLetter[i][j]) == 0)
			{
				copyNFA(NFA, copy, NFA[state].statesWithLetter[i][j], states, copySize);
			}
		}
	}
}

vector<NFAState> getNFA(int start, vector<NFAState> &NFA, int &NFASize)
{
	vector<NFAState> copy;
	for (int i = 0; i < NFA.size(); i++) {
		copy.push_back(q);
	}
	set<int> states;
	int copySize = 0;
	copyNFA(NFA, copy, start, states, copySize);
	for (int i = 0; i < copySize; ++i)
	{
		NFASize--;
		NFA.pop_back();
	}
	return copy;
}

void pushDFA(vector<DFAState> DFA, vector<NFAState> &NFA, int &NFASize, int finalState)
{
	for (int i = 0; i < DFA.size(); ++i)
	{
		NFA.push_back(q);
		NFA[NFASize].isFinal = false;
		if (!DFA[i].isFinal)
		{
			NFA[NFASize].statesWithEpsilon.push_back(finalState);
		}

		for (int j = 0; j < alphabetSize; ++j)
		{
			NFA[NFASize].statesWithLetter[j].push_back(DFA[i].stateWithLetter[j] + finalState + 1);
		}
		NFASize++;
	}
}

void complement(vector<NFAState> &NFA, stack<int> &endpoints, int &NFASize)
{
	int q2 = endpoints.top();
	endpoints.pop();
	int q1 = endpoints.top();
	endpoints.pop();

	vector<NFAState> negation = getNFA(q1, NFA, NFASize);

	negation[q2].isFinal = true;

	vector<DFAState> detNegation = determinize(negation, q1);
	NFA.push_back(q);
	NFASize++;
	int s = NFASize;
	endpoints.push(NFASize);
	endpoints.push(NFASize - 1);
	pushDFA(detNegation, NFA, NFASize, NFASize - 1);
}

void intersect(vector<NFAState> &NFA, int &NFASize, stack<int> &endpoints)
{
	int q4 = endpoints.top(); endpoints.pop();
	int q3 = endpoints.top(); endpoints.pop();
	int q2 = endpoints.top(); endpoints.pop();
	int q1 = endpoints.top(); endpoints.pop();

	vector<NFAState> B = getNFA(q3, NFA, NFASize);
	B[q4].isFinal = true;
	vector<DFAState> detB = determinize(B, q3);

	vector<NFAState> A = getNFA(q1, NFA, NFASize);
	A[q2].isFinal = true;
	vector<DFAState> detA = determinize(A, q1);

	NFA.push_back(q);
	NFASize++;
	int s1 = NFASize;
	pushDFA(detA, NFA, NFASize, NFASize - 1);
	NFA.push_back(q);
	endpoints.push(NFASize);
	NFA[NFASize].statesWithEpsilon.push_back(s1);
	endpoints.push(s1 - 1);
	NFASize++;

	NFA.push_back(q);
	NFASize++;
	int s2 = NFASize;
	pushDFA(detB, NFA, NFASize, NFASize - 1);
	NFA.push_back(q);
	endpoints.push(NFASize);
	NFA[NFASize].statesWithEpsilon.push_back(s2);
	endpoints.push(s2 - 1);
	NFASize++;

	unite(NFA, NFASize, endpoints);

	complement(NFA, endpoints, NFASize);
}

void difference(vector<NFAState> &NFA, int &NFASize, stack<int> &endpoints)
{
	int q4 = endpoints.top(); endpoints.pop();
	int q3 = endpoints.top(); endpoints.pop();

	vector<NFAState> B = getNFA(q3, NFA, NFASize);
	B[q4].isFinal = true;
	vector<DFAState> detB = determinize(B, q3);

	NFA.push_back(q);
	NFASize++;
	int s2 = NFASize;
	pushDFA(detB, NFA, NFASize, NFASize - 1);
	NFA.push_back(q);
	endpoints.push(NFASize);
	NFA[NFASize].statesWithEpsilon.push_back(s2);
	endpoints.push(s2 - 1);
	NFASize++;

	intersect(NFA, NFASize, endpoints);
}

void changeTransitions(int q1, vector<NFAState> NFA, vector<NFAState> &result, set<int> &visited)
{
	visited.insert(q1);
	for (int i = 0; i < NFA[q1].statesWithEpsilon.size(); ++i)
	{
		result[NFA[q1].statesWithEpsilon[i]].statesWithEpsilon.push_back(q1);
		if (visited.count(NFA[q1].statesWithEpsilon[i]) == 0)
		{
			changeTransitions(NFA[q1].statesWithEpsilon[i], NFA, result, visited);
		}
	}
	for (int i = 0; i < alphabetSize; ++i)
	{
		for (int j = 0; j < NFA[q1].statesWithLetter[i].size(); ++j)
		{
			result[NFA[q1].statesWithLetter[i][j]].statesWithLetter[i].push_back(q1);
			if (visited.count(NFA[q1].statesWithLetter[i][j]) == 0)
			{
				changeTransitions(NFA[q1].statesWithLetter[i][j], NFA, result, visited);
			}
		}
	}
}

void copyToNFA(int q2, vector<NFAState> &NFA, vector<NFAState> result, set<int> &visited)
{
	visited.insert(q2);
	for (int i = 0; i < result[q2].statesWithEpsilon.size(); ++i)
	{

		if (visited.count(result[q2].statesWithEpsilon[i]) == 0)
		{
			copyToNFA(result[q2].statesWithEpsilon[i], NFA, result, visited);
		}
	}
	for (int i = 0; i < alphabetSize; ++i)
	{
		for (int j = 0; j < result[q2].statesWithLetter[i].size(); ++j)
		{
			if (visited.count(result[q2].statesWithLetter[i][j]) == 0)
			{
				copyToNFA(result[q2].statesWithLetter[i][j], NFA, result, visited);
			}
		}
		NFA[q2].statesWithLetter[i] = result[q2].statesWithLetter[i];
	}
	NFA[q2].statesWithEpsilon = result[q2].statesWithEpsilon;

}

void reverse(vector<NFAState> &NFA, stack<int> &endpoints, int NFASize)
{
	int q2 = endpoints.top();
	endpoints.pop();
	int q1 = endpoints.top();
	endpoints.pop();
	vector<NFAState> result;
	set<int> visited;
	for (int i = 0; i < NFASize; ++i)
	{
		result.push_back(q);
	}
	changeTransitions(q1, NFA, result, visited);
	set<int> visited3;
	copyToNFA(q2, NFA, result, visited3);
	endpoints.push(q2);
	endpoints.push(q1);
}

pair<vector<NFAState>, int> constructNFA(string expression)
{
	vector<NFAState> NFA;
	int NFASize = 0;
	stack<int> endpoints;
	for (unsigned int i = 0; i < expression.size(); ++i)
	{
		switch (expression[i])
		{
		case '*': kleeneStar(NFA, NFASize, endpoints); break;
		case '.': concatenate(NFA, endpoints); break;
		case '+': unite(NFA, NFASize, endpoints); break;
		case '-': difference(NFA, NFASize, endpoints); break;
		case '&': intersect(NFA, NFASize, endpoints); break;
		case '^': complement(NFA, endpoints, NFASize); break;
		case '<': reverse(NFA, endpoints, NFASize); break;
		case '0': letter(26, NFA, NFASize, endpoints);
		case '1': letter(27, NFA, NFASize, endpoints);
		case '2': letter(28, NFA, NFASize, endpoints);
		case '3': letter(29, NFA, NFASize, endpoints);
		case '4': letter(30, NFA, NFASize, endpoints);
		case '5': letter(31, NFA, NFASize, endpoints);
		case '6': letter(32, NFA, NFASize, endpoints);
		case '7': letter(33, NFA, NFASize, endpoints);
		case '8': letter(34, NFA, NFASize, endpoints);
		case '9': letter(35, NFA, NFASize, endpoints);
		case ' ': letter(36, NFA, NFASize, endpoints);
		case ',': letter(37, NFA, NFASize, endpoints);
		case '$': letter(38, NFA, NFASize, endpoints);
		default: letter(expression[i] - 'A', NFA, NFASize, endpoints);
		}
	}
	int f; //final state of NFA
	int q0; //initial state of NFA
	f = endpoints.top(); endpoints.pop();
	q0 = endpoints.top(); endpoints.pop();
	NFA[f].isFinal = true;
	return make_pair(NFA, q0);
}

pair<int, vector<DFAState> > minimize(vector<DFAState> DFA)
{
	vector<int> groupNumber(DFA.size());
	vector<vector<int> > partitions(2, vector<int>());

	partitions[0].push_back(0);
	for (int i = 1; i < (int)groupNumber.size(); ++i)
	{
		if (DFA[i].isFinal == DFA[0].isFinal)
		{
			groupNumber[i] = 0;
			partitions[0].push_back(i);
		}
		else
		{
			groupNumber[i] = 1;
			partitions[1].push_back(i);
		}
	}

	if (!partitions[1].size())
	{
		partitions.erase(partitions.end());
	}
	bool canPart = true;
	int initialState = 0;
	while (canPart)
	{
		canPart = false;

		for (int i = 0; i < partitions.size(); ++i)
		{
			for (int j = 0; j < alphabetSize; ++j)
			{
				vector<pair<int, int> > groupState(partitions[i].size());
				for (int k = 0; k<partitions[i].size(); ++k)
				{
					if (DFA[partitions[i][k]].stateWithLetter[j] >= 0)
					{
						groupState[k] = make_pair(groupNumber[DFA[partitions[i][k]].stateWithLetter[j]], partitions[i][k]);
					}
					else
					{
						groupState[k] = make_pair(-1, partitions[i][k]);
					}
				}
				sort(groupState.begin(), groupState.end());

				if (groupState[0].first != groupState[groupState.size() - 1].first)
				{
					canPart = true;

					int k, m = partitions.size() - 1;

					partitions[i].clear();
					partitions[i].push_back(groupState[0].second);
					for (k = 1; k<groupState.size() && (groupState[k].first == groupState[k - 1].first); k++)
					{
						partitions[i].push_back(groupState[k].second);
					}

					while (k<groupState.size())
					{
						if (groupState[k].first != groupState[k - 1].first)
						{
							partitions.push_back(vector<int>());
							m++;
						}
						groupNumber[groupState[k].second] = m;
						partitions[m].push_back(groupState[k].second);
						k++;
					}
				}
			}
		}
	}

	for (int i = 0; i < partitions.size(); ++i)
	{
		for (int j = 0; j < partitions[i].size(); ++j)
		{
			if (partitions[i][j] == 0)
			{
				initialState = i;
				break;
			}
		}
	}

	vector<DFAState> minDFA(partitions.size());

	for (int i = 0; i<(int)partitions.size(); ++i)
	{
		for (int j = 0; j < alphabetSize; ++j)
		{
			if (DFA[partitions[i][0]].stateWithLetter[j] >= 0)
			{
				minDFA[i].stateWithLetter[j] = groupNumber[DFA[partitions[i][0]].stateWithLetter[j]];
			}
			else
			{
				minDFA[i].stateWithLetter[j] = -1;
			}
		}

		minDFA[i].isFinal = DFA[partitions[i][0]].isFinal;
	}

	return make_pair(initialState, minDFA);
}



int main()
{
	string expression;
	string month29 = "FE.B.R.U.A.R.Y.";
	string month30 = "AP.R.I.L.JU.N.E.SE.P.T.E.M.B.E.R.NO.V.E.M.B.E.R.+++";
	string month31 = "JANUARY......MARCH....MAY..JULY...AUGUST.....OCTOBER......DECEMBER.......++++++";
	string month = month29 + month30 + month31 + "++";
	string oneToNine = "123456789++++++++";
	string zeroToNine = "0" + oneToNine + "+";
	string date = oneToNine + "12+" + zeroToNine + ".30.31.+++";
	string year = oneToNine + "A*AA*.-" + zeroToNine + "+" + "A*AA*.-" + zeroToNine + "+" + "A*AA*.-" + zeroToNine + "+...";
	string allDates = month + " ." + date + ", ..." + year + "$..";
	string sigma = "AB+C+D+E+F+G+H+I+J+K+L+M+N+O+P+Q+R+S+T+U+V+W+X+Y+Z+0+1+2+3+4+5+6+7+8+9+ +,+$+";
	string maxDays30 = sigma + "*" + sigma + "*" + month29 + " 30...." + sigma + "*.-";
	string maxDays31 = sigma + "*" + sigma + "*" + month30 + " 31...." + sigma + "*.-";
	string maxDays = maxDays30 + maxDays31 + "&";
	string maxDates = allDates + maxDays + "&";
	string even = "02468++++";
	string odd = "13579++++";
	string N = even + odd + "+";
	string div4 = "48+" + N + "*" + even + "048++." + odd + "26+.+.+";
	string leapYear = div4 + N + N + "*." + div4 + "-" + "00..-";
	string leapDates = sigma + "*" + month29 + " 29, ......" + leapYear + "$" + sigma + "*..^.^";
	string validDates = allDates + maxDays + leapDates + "&&";
	/*cout << "Regular expression in reverse polish notation: ";
	cin >> expression;*/

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	pair<vector<NFAState>, int> NFAconstruction = constructNFA(month29);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	auto duration = duration_cast<microseconds>(t2 - t1).count();
	cout << duration << endl;

	int q0 = NFAconstruction.second;
	vector<NFAState> NFA = NFAconstruction.first;

	high_resolution_clock::time_point t3 = high_resolution_clock::now();
	vector<DFAState> DFA = determinize(NFA, q0);
	high_resolution_clock::time_point t4 = high_resolution_clock::now();

	auto duration2 = duration_cast<microseconds>(t4 - t3).count();
	cout << duration2 << endl;

	bool hasFinals = false;
	for (int i = 0; i < DFA.size(); ++i)
	{
		if (DFA[i].isFinal)
		{
			hasFinals = true;
			break;
		}
	}

	if (hasFinals)
	{
		high_resolution_clock::time_point t5 = high_resolution_clock::now();
		pair<int, vector<DFAState> > result = minimize(DFA);
		high_resolution_clock::time_point t6 = high_resolution_clock::now();

		auto duration3 = duration_cast<microseconds>(t6 - t5).count();
		cout << duration3 << endl;

		vector<DFAState>  minDFA = result.second;
		int initialState = result.first;
		cout << "Initial state: " << initialState << endl;

		high_resolution_clock::time_point t7 = high_resolution_clock::now();
		print(minDFA);
		high_resolution_clock::time_point t8 = high_resolution_clock::now();

		auto duration4 = duration_cast<microseconds>(t8 - t7).count();
		cout << duration4 << endl;
	}
	else
	{
		cout << "The automaton doesn't accept anything." << endl;
	}

	return 0;
}
