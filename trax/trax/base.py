

CHECKPOINT_FILE_MODE = 'wb'
LOG_FILE_MODE = 'wb'

class AbstractTransactional(object):

	def __init__(self, checkpoint='transactional.cpt', log='transactional.log', checkpoint_mode=CHECKPOINT_FILE_MODE, log_mode=LOG_FILE_MODE):
		assert type(checkpoint) is str, type(checkpoint)
		assert type(log)        is str, type(log)

		self._checkpoint_path = checkpoint
		self._log_path        = log
		self._cpt_mode        = checkpoint_mode
		self._log_mode        = log_mode
		self._checkpoint_fd   = None
		self._log_fd          = None

	@property
	def cpt_path(self): return self._checkpoint_path

	@property
	def log_path(self): return self._log_path

	def close(self):
		self._cpt_close()
		self._log_close()

	def __enter__(self):
		return self

	def __exit__(self, exc_type, exc_value, traceback):
		self.close()

	def _log_close(self):
		if self._log_fd is not None:
			self._log_fd.close()
			self._log_fd = None

	def _cpt_close(self):
		if self._checkpoint_fd is not None:
			self._checkpoint_fd.close()
			self._checkpoint_fd = None

	def _log_open(self):
		if self._log_fd is None:
			self._log_fd = open(self.log_path, self._log_mode)

	def _cpt_open(self):
		if self._checkpoint_fd is None:
			self._checkpoint_fd = open(self.cpt_path, self._cpt_mode)

	def checkpoint(self, value):
		self._log_close()
		self._cpt_open()
		self._impl_checkpoint(self._checkpoint_fd, value)
		self._cpt_close()

	def _impl_checkpoint(self, fd, value):
		raise NotImplementedError

	def log(self, value):
		self._log_open()
		self._impl_log(self._log_fd, value)

	def _impl_log(self, fd, value):
		raise NotImplementedError

	def recover(self, checkpoint_handler=None, log_handler=None):
		"""
		checkpoint_handler :: FileHandle -> IO a
		log_handler        :: a -> FileHandle -> IO a
		"""
		self._cpt_close()
		self._log_close()
		with self._impl_cpt_recover_open() as fd:
			obj = checkpoint_handler(fd)
		with self._impl_log_recover_open() as fd:
			obj = log_handler(obj, fd)
		return obj

	def _impl_cpt_recover_open(self):
		"""
		:: IO FileHandle
		"""
		raise NotImplementedError

	def _impl_log_recover_open(self):
		"""
		:: IO FileHandle
		"""
		raise NotImplementedError
