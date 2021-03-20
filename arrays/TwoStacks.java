package arrays;

public class TwoStacks {
	
	int arr[];
	int top1 = -1;
	int top2 ;
	
	TwoStacks(int size)
	{
		arr = new int[size];		
		top2 = size;
	}
	
	void push1(int data)
	{
		if(top2 == top1+1)
			System.out.println("Stack is overflowed");
		
		else
		{
			arr[++top1] = data;
		}
	}
	
	void push2(int data)
	{
		if(top2 == top1+1)
			System.out.println("Stack is overflowed");
		
		else
		{
			arr[--top2] = data;
		}
	}
	
	void pop1()
	{
		if(top1 == -1)
			System.out.println("stack is underflowed");
		
		else
		{
			top1--;
		}
	}
	
	void pop2()
	{
		if(top2 == arr.length)
			System.out.println("stack is underflowed");
		
		else
		{
			top2--;
		}
	}
	
	void traverse1()
	{
		if(top1 == -1)
			System.out.println("Stack is empty");
		else
		{
			for(int i=top1;i>=0;i--)
			{
				System.out.print(arr[i]+" ");
			}
			System.out.println();
		}
	}
	
	void traverse2()
	{
		if(top2 == arr.length)
			System.out.println("Stack is empty");
		else
		{
			for(int i=top2;i<arr.length;i++)
			{
				System.out.print(arr[i]+" ");
			}
			System.out.println();
		}
	}
	
	public static void main(String[] args) {
		
		TwoStacks st =new TwoStacks(8);
		
		st.push1(10);
		st.push1(20);
		st.push1(30);
		st.push1(40);
		st.push2(110);
		st.push2(120);
		st.push2(130);
		st.push2(140);
		st.push2(150);
		
		st.traverse1();
		st.traverse2();
		
		
		

	}

}
