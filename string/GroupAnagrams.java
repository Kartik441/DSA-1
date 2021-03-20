import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;


public class GroupAnagrams {
	
	// Other approaches on gfg

	// efficient approach using trie Advanced Data Structure
	
	// Hashmap solution  TC: O(nm)  n words and each word has maximum m characters
						//SC: O(n+m)
	static ArrayList<ArrayList<String>> solve(ArrayList<String> list)
	{
//		HashMap<HashMap<Character, Integer>, ArrayList<String>> map = new HashMap<HashMap<Character, Integer>,
//				ArrayList<String>>();
		HashMap<HashMap<Character, Integer>, ArrayList<String>> map = new HashMap<>();
		
		for(String str : list)
		{
			HashMap<Character, Integer> tempMap = new HashMap<>();
			// counting frequency of characters present in a string
			for(int i = 0;i < str.length() ; i++)
			{
				if(tempMap.containsKey(str.charAt(i)))
				{
					int x = tempMap.get(str.charAt(i));
					tempMap.put(str.charAt(i), ++x);
				}
				
				else
				{
					tempMap.put(str.charAt(i), 1);
				}
			}
			
			if(map.containsKey(tempMap))
				map.get(tempMap).add(str);
			else
			{
				ArrayList<String> tempList = new ArrayList<String>();
				tempList.add(str);
				map.put(tempMap, tempList);
			}
		}
		ArrayList<ArrayList<String>> res = new ArrayList<>();
		for(HashMap<Character, Integer> temp : map.keySet())
			res.add(map.get(temp));
		return res;
	}
	////////////////////////////////////////////////////////////////////////////////////////////
	// Using Trie
	// Insert the sorted order of each word in a trie
	// Since all the anagrams will end at the same leaf node, we can start a linked list at the leaf
	// nodes where each node represents the index of the original array of words
	static class TNode
	{
		boolean isEnd;
		TNode child[] = new TNode[26];
		LinkedList<Integer> list;
		
		TNode()
		{
			isEnd =false;
			list = new LinkedList<>();
			
		}
		
	}
	
	static void insert(TNode root, String s, int orInd)
	{
		int n = s.length();
		TNode temp = root;
		for(int i=0;i<n;i++)
		{
			int index = s.charAt(i) - 'a';
			
			if(temp.child[index] == null)
				temp.child[index] = new TNode();
			
			temp = temp.child[index];
		}
		
		if(temp.isEnd == true)
		{
			temp.list.add(orInd);
		}
		else
		{
			temp.isEnd = true;
			temp.list.add(orInd);
		}
		
			
	}
	
	static void printAnagrams(ArrayList<String> arr)
	{
		int n = arr.size();
		TNode root = new TNode();
		for(int i=0;i<n;i++)
		{
			char[] buffer = arr.get(i).toCharArray();
			// Sort the buffer
			Arrays.sort(buffer);
			insert(root, new String(buffer), i);
		}
		
		// traverse the trie and print the answer
		printAnagrams(root, arr);
	}
	
	static void printAnagrams(TNode root, ArrayList<String> arr)
	{
		if(root == null)
			return;
		
		if(root.isEnd)
		{
			// traverse the list
			for(int index : root.list)
				System.out.print(arr.get(index)+" ");
			System.out.println();
		}
		
		// traverse all the child nodes
		for(int i=0;i<26;i++)
			printAnagrams(root.child[i], arr);
	}

	public static void main(String[] args) {
		ArrayList<String> list = new ArrayList<>();
		list.add("cat");
		list.add("dog");
		list.add("ogd");
		list.add("god");
		list.add("atc");
		
		System.out.println(solve(list));
		System.out.println();
		printAnagrams(list);

	}

}
