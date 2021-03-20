package ads;

// can be used to print all unique rows in a given binary matrix
// the trie node will contain child array of 2 size, one for 0 and other for 1 in matrix
// insert each row in the trie, if it is already present skip it otherwise print it

class Trnode
{
	char data;
	int we;                // word end
	Trnode child[];
	Trnode()
	{
		this.we = -1;
		child = new Trnode[26];        // For 26 characters
	}
	Trnode(char data, int we)
	{
		this.data = data;
		this.we = we;
		child = new Trnode[26];        // For 26 characters
	}
}

public class Trie {
	
	Trnode root;
	
	static void insert(Trnode root, String s)
	{
		int n = s.length();
		Trnode temp = root;
		
		for(int i=0;i<n;i++)
		{
			char c = s.charAt(i);
			int index = c - 'a';
			if(temp.child[index] == null)
			{
				temp.child[index] = new Trnode(c, 0);
			}
			
			temp = temp.child[index];
			
		}
		temp.we++;
	}
	
	
	static boolean search(Trnode root, String s)
	{
		int n = s.length();
		Trnode temp = root;
		for(int i=0;i<n;i++)
		{
			char c = s.charAt(i);
			int index = c - 'a';
			
			if(temp.child[index] == null)
				return false;
			
			temp = temp.child[index];
		}
		return temp != null && temp.we != 0; 
 
		
	}
	
	static void delete(Trnode root, String s)
	{
		int n = s.length();
		Trnode temp = root;
		for(int i=0;i<n;i++)
		{
			char c = s.charAt(i);
			int index = c - 'a';
			
			if(temp.child[index] == null)
			{
				System.out.println("String not present");
				return;
			}
			
			temp = temp.child[index];
		}
		if(temp != null)
			temp.we--;
		else
			System.out.println("String not present");
		
	}
	

	public static void main(String[] args) {
		Trie t = new Trie();
		t.root = new Trnode();
		
		
		String keys[] = {"the", "a", "there", "answer", "any", "by", "bye", "their"};
		for(int i=0;i<keys.length;i++)
		{
			insert(t.root, keys[i]);
		}
		
		System.out.println(search(t.root, "the"));
		System.out.println(search(t.root, "these"));
		System.out.println(search(t.root, "their"));
		System.out.println(search(t.root, "thaw"));
		
		delete(t.root, "their");
		System.out.println(search(t.root, "their"));
		
	}

}
