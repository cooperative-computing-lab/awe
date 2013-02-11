
import trax

import os, random, unittest


import cPickle as pickle



class  TestSimple(unittest.TestCase):

	def clean(self):
		for p in ['transactional.cpt', 'transactional.log']:
			if os.path.exists(p):
				os.unlink(p)

	def setUp(self):
		self.clean()

	def tearDown(self):
		self.clean()


	def test_example(self):

		trx = trax.SimpleTransactional(pickleprotocol=-1)
		state = []
		for i in xrange(1000):
			v = random.randint(0,100000)
			state.append(v)
			trx.log(v)
			if i % 50 == 0:
				trx.checkpoint(state)

		recovered = trx.recover(lambda o, v: o.append(v))

		self.assertTrue( state == recovered )


if __name__ == '__main__':
	TestExample().test_example()
