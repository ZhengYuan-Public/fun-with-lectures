



class ClassName(object):
	
	ClassVar = var_1

	def __init__(self, arg1, arg2, ...):
		self.arg1 = arg1
		self.arg2 = arg2
		...

	def instance_method_name(self):
		func(self.arg1, self.arg2, self.ClassVar, ...)

	@classmethod # Can also be uses as Alternative Constructor
	def class_method_name(cls, arg1, arg2, ...):
		func(cls.arg1, cls,arg2, ...)
 
	@staticmethod
	def static_method_name(arg_from_outside_class):
		func(arg_from_outside_class)

class InheritedSubClass(ClassName):
	def __init__(self, arg1, arg2, ..., new_arg1, new_arg2, ...):
		super().__init__(arg1, arg2, ...)
		# Or use following
		ClassName.__init__(self, arg1, arg2, ...)
		self.new_arg1  = new_arg1
		self.new_arg2 = new_arg2
		...

	def child_class_method_name(self):
		func(self.new_arg1, self.new_arg2)







class Name:
	def __init__(self, first, last):
		self.first = first
		self.last = last

	def full_name(self):
		return '{} {}'.format(self.first, self.last)

name_1 = Name('Zheng', 'Yuan')
name_1.full_name()




class Name:
	def __init__(self, first, last):
		self.first = first
		self.last = last
	
	@property
	def fullname(self):
		return '{} {}'.format(self.first, self.last)

	@fullname.setter
	def fullname(self, name):
		first, last = name.split(' ')
		self.first = first
		self.last = last

	@fullname.deleter
	def fullname(self)
		print('Delete Fullname!')
		self.first = None
		self.last = None

name_1 = Name('Zheng', 'Yuan')

name_1.fullname = 'John Doe'

name_1.full_name

del name_1.fullname




	



		




